(ns patternq.paging-test
  "Set paging and counting against a fake /datoms: with cursor paging and
  in-place counts (the 2026-11 service) or offsets only (the one before)."
  (:require [clojure.test :refer [deftest is testing]]
            [patternq.db :as pdb]
            [patternq.dataset :as pqd]))

(defn- fake-set [n cursor? countable?]
  (let [ids (vec (range 1000 (+ 1000 n)))
        calls (atom 0)]
    {:ids ids
     :calls calls
     :datoms (fn [_db _index _components {:keys [offset limit after] :or {offset 0}}]
               (assert (<= offset 100000) "offset past the service cap")
               (swap! calls inc)
               (let [start (if after (parse-long after) offset)
                     page (subvec ids (min start n) (min (+ start limit) n))]
                 (with-meta (mapv (fn [v] {:e 1 :v v}) page)
                            {:patternq/cursor? cursor?
                             :patternq/next (when (and cursor? (= (count page) limit))
                                              (str (+ start (count page))))})))
     :datoms-count (fn [_db _index _components cap]
                     (when countable? [(min n (inc cap)) (<= n cap)]))}))

(deftest set-paging-with-and-without-cursor
  (doseq [cursor? [true false] n [0 7 10000 25001 150000]]
    (testing (str "cursor? " cursor? ", n " n)
      (let [{:keys [ids datoms datoms-count]} (fake-set n cursor? false)]
        (with-redefs [pdb/datoms datoms pdb/datoms-count datoms-count]
          (is (= [(vec (take 20000 ids)) (<= n 20000)]
                 (#'pqd/set-measurement-ids :db 1 20000)))
          (is (= (if (<= n 100000) [n true] [100000 false])
                 (#'pqd/count-set-measurements :db 1 100000)
                 ;; without in-place counts, a bigger cap stops at the 100000 offset cap
                 (#'pqd/count-set-measurements :db 1 1000000))))))))

(deftest in-place-count-skips-paging
  (let [{:keys [calls datoms datoms-count]} (fake-set 150000 true true)]
    (with-redefs [pdb/datoms datoms pdb/datoms-count datoms-count]
      (is (= [100000 false] (#'pqd/count-set-measurements :db 1 100000)))
      (is (zero? @calls)))))

(deftest cursor-page-from-an-old-instance-is-reread-by-offset
  (let [ids (vec (range 1000 26001))]
    (with-redefs [pdb/datoms (fn [_ _ _ {:keys [offset limit after] :or {offset 0}}]
                               (let [new? (and (nil? after) (zero? offset))
                                     start (if after 0 offset)
                                     page (subvec ids (min start 25001) (min (+ start limit) 25001))]
                                 (with-meta (mapv (fn [v] {:v v}) page)
                                            {:patternq/cursor? new?
                                             :patternq/next (when new? (str (+ start (count page))))})))]
      (is (= [ids true] (#'pqd/set-measurement-ids :db 1 30000))))))
