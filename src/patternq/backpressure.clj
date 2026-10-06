(ns patternq.backpressure
  "Client side of the commons API back-pressure contract (unify-central
  dev-docs/API_BACKPRESSURE.md), the same in the R, Python, Clojure and Julia
  libraries.

  - 429 and 503 are retried, up to :max-retries times, when the body is a
    retryable application/problem+json, or is not problem+json at all (e.g.
    from a load balancer). Every other status, 502 and 504 included, is
    returned as is.
  - The wait is Retry-After when the server sends it, otherwise exponential
    backoff with full jitter: random(0, min(:max-backoff, 2^attempt)) seconds.
  - Self-throttling: a RateLimit header with r=0 makes the next call wait t
    seconds, and at most :max-concurrency /query and /datoms calls run at once
    in this process (the server's per-key cap is 4).
  - Retries are logged at INFO on the \"patternq\" System.Logger, with the
    problem type and the wait.

  Presigned S3 downloads are not API calls and don't go through here."
  (:require [charred.api :as json]
            [clojure.string :as str])
  (:import (java.lang System$Logger System$Logger$Level)
           (java.net.http HttpResponse)
           (java.util.concurrent Semaphore)))

(set! *warn-on-reflection* true)

(def problem-prefix "urn:pattern-data-commons:problem:")
(def default-max-timeout-ms 120000)

(def ^:private default-policy {:max-retries 5 :max-backoff 60 :max-concurrency 4})
(defonce ^:private policy (atom default-policy))
(defonce ^:private slots (atom (Semaphore. 4 true)))
(defonce ^:private next-allowed-ms (atom 0))
(defonce ^:private max-timeout (atom default-max-timeout-ms))

(def ^:private ^System$Logger logger (System/getLogger "patternq"))

(defn log!
  [level ^String msg]
  (.log logger (case level :info System$Logger$Level/INFO :warning System$Logger$Level/WARNING) msg))

(defn- sleep! [ms] (Thread/sleep (long ms)))

(defn retry-policy
  "The current back-pressure policy (see set-retry-policy!)."
  []
  @policy)

(defn set-retry-policy!
  "Set how this process handles API back pressure. Keys left out keep their
  current value. Defaults: :max-retries 5 retries per call, :max-backoff 60
  seconds (cap on the jittered backoff when the server sends no Retry-After),
  :max-concurrency 4 concurrent /query and /datoms calls (the server's per-key
  cap). Returns the previous policy."
  [{:keys [max-retries max-backoff max-concurrency] :as opts}]
  (when (and max-concurrency (< max-concurrency 1))
    (throw (ex-info "max-concurrency must be at least 1" opts)))
  (let [[old new] (swap-vals! policy merge (select-keys opts [:max-retries :max-backoff :max-concurrency]))]
    (when (not= (:max-concurrency old) (:max-concurrency new))
      (reset! slots (Semaphore. (int (:max-concurrency new)) true)))
    old))

(defn max-timeout-ms
  "The server's cap on requested query timeouts, from the last
  PDC-Query-Max-Timeout-Ms header seen (120000 until one is)."
  []
  @max-timeout)

(defn- header [^HttpResponse resp name]
  (.orElse (.firstValue (.headers resp) name) nil))

(defn body-string
  "A response body (String or byte[]) as a string."
  [^HttpResponse resp]
  (let [b (.body resp)]
    (if (bytes? b) (String. ^bytes b "UTF-8") (str b))))

(defn- problem [^HttpResponse resp]
  (when (str/starts-with? (or (header resp "Content-Type") "") "application/problem+json")
    (let [p (try (json/read-json (body-string resp)) (catch Exception _ nil))]
      (when (map? p) p))))

(defn- retry-after-s [resp problem]
  (or (some-> (header resp "Retry-After") str/trim parse-double (max 0.0))
      (when (number? (get problem "retry_after")) (max 0.0 (double (get problem "retry_after"))))))

(defn- note-headers! [resp]
  (when-let [rl (header resp "RateLimit")]
    (let [r (some->> (re-find #"\br=(\d+)" rl) second parse-long)
          t (some->> (re-find #"\bt=(\d+)" rl) second parse-long)]
      (when (and (= 0 r) t)
        (swap! next-allowed-ms max (+ (System/currentTimeMillis) (* 1000 t))))))
  (when-let [mt (some-> (header resp "PDC-Query-Max-Timeout-Ms") str/trim parse-long)]
    (reset! max-timeout mt)))

(defn- wait-for-rate-limit! []
  (let [wait (- @next-allowed-ms (System/currentTimeMillis))]
    (when (pos? wait)
      (log! :info (format "patternq: rate limit reached (RateLimit r=0); waiting %.1f s before the next call"
                          (/ wait 1000.0)))
      (sleep! wait))))

(defn send!
  "Send an API request (do-request performs it and returns the HttpResponse)
  under the back-pressure contract. limited?: the call counts toward the
  per-key concurrency cap (/query, /datoms). Returns the final response,
  which may still be an error for the caller to throw; throws ex-info with
  :patternq/throttled true when throttling outlasts :max-retries."
  ^HttpResponse [do-request what limited?]
  (loop [attempt 0]
    (wait-for-rate-limit!)
    (let [^HttpResponse resp (if limited?
                               (let [^Semaphore s @slots]
                                 (.acquire s)
                                 (try (do-request) (finally (.release s))))
                               (do-request))
          _ (note-headers! resp)
          status (.statusCode resp)
          p (when (#{429 503} status) (problem resp))]
      (if (or (not (#{429 503} status))
              (and p (not (true? (get p "retryable")))))
        resp
        (let [ptype (get p "type")
              label (or ptype (str "HTTP " status))
              {:keys [max-retries max-backoff]} @policy]
          (when (>= attempt max-retries)
            (throw (ex-info (format "%s throttled by the commons API (%s) and still throttled after %d retries: %s"
                                    what label attempt
                                    (or (get p "detail") (let [b (str/trim (body-string resp))]
                                                           (subs b 0 (min 500 (count b))))))
                            {:patternq/throttled true :problem-type ptype :status status})))
          (let [attempt (inc attempt)
                wait-s (or (retry-after-s resp p)
                           (* (rand) (min (double max-backoff) (Math/pow 2 attempt))))]
            (log! :info (format "patternq: %s throttled (%s, HTTP %d); retry %d of %d in %.1f s"
                                what label status attempt max-retries wait-s))
            (sleep! (Math/round (* 1000.0 (double wait-s))))
            (recur attempt)))))))
