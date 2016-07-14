(ns clj-reciprocal-blast.core
  (:require [clj-fasta.core :as fa]
            [clj-blast.core :as bl]
            [biodb.core :as bdb]
            [me.raynes.fs :as fs]
            [clojure.java.io :as io]
            [clojure.string :as st]))

(defn- make-wd
  [wd]
  (if-not wd
    (loop [dirname (str "./" "rec-blast-" (java.util.UUID/randomUUID) "/")
           c (atom 0)]
      (cond (not (fs/directory? dirname))
            (if (fs/mkdir dirname)
              dirname)
            (> @c 10)
            (throw (Exception. "Cannot create working directory."))
            :else
            (recur (str wd "/" "rec-blast-" (str (java.util.UUID/randomUUID))) (swap! c inc))))
    wd))

(defn- blasts-match?
  [iteration db relax]
  (let [seqid (bl/accession iteration)
        dbhit (->> (bdb/get-sequences db :sequences :blast
                                      [(-> (bl/hit-seq iteration) first :Hit_id)])
                   first
                   bl/hit-seq
                   (take relax)
                   (map :Hit_id))]
    (some #(= seqid %) dbhit)))

(defn- first-db-type
  [program]
  (condp = program
    "blastp" "prot"
    "blastx" "prot"
    "blastn" "nucl"
    "tblastn" "nucl"))

(defn- second-db-type
  [program]
  (condp = program
    "blastp" "prot"
    "blastx" "nucl"
    "blastn" "nucl"
    "tblastn" "prot"))

(defn- second-program
  [program]
  (condp = program
    "blastx" "tblastn"
    "blastp" "blastp"
    "blastn" "blastn"
    "tblastn" "blastx"))

(defn- significant?
  [i evalue query-coverage]
  (if-let [h (-> (bl/hit-seq i :evalue evalue) first)]
    (if (>= (/ (bl/query-length i) (-> (:hsps h) first :Hsp_align-len))
            query-coverage)
      (:Hit_id h))))

(defn- sign-hits->file
  [{:keys [blast1 wdir evalue database program query-coverage] :as m}]
  (let [of (str (fs/file wdir "bl1-significant-hits.fasta"))]
    (-> (mapcat
         #(with-open [r (io/reader %)]
            (doall
             (->> (bl/iteration-seq r)
                  (map (fn [x] (significant? x evalue query-coverage)))
                  (remove nil?))))
         blast1)
        set
        (bl/blastdb->file database of (first-db-type program)))
    (assoc m :sign-fasta of)))

(defn- first-blasts
  [{:keys [sfile program database wdir bparams] :as m}]
  (let [res (->> (assoc m :blast1
                        (bl/blast-file sfile program database (str (fs/file wdir "blast1"))
                                       :params (merge bparams {"-max_target_seqs" "1"})))
                 sign-hits->file)]
    res))

(defn- second-blasts
  [{:keys [blast1 program sfile wdir bparams sign-fasta] :as m}]
  (let [res (assoc m :blast2
                   (bl/blast-file sign-fasta
                                  (second-program program)
                                  (bl/create-blastdb-file sfile (second-db-type program))
                                  (str (fs/file wdir "blast2"))
                                  :params bparams))]
    res))

(defn- save-blasts
  [{:keys [blast2] :as m}]
  (let [db (bdb/db-spec {:dbname (str (fs/file
                                       (fs/parent (first blast2))
                                       (first (st/split (fs/name (first blast2)) #"-")))
                                      ".sqlite")
                         :dbtype :sqlite})]
    (bdb/create-table! db :sequences :blast)
    (doseq [f blast2]
      (with-open [r (io/reader f)]
        (bdb/insert-sequences! db :sequences :blast (bl/iteration-seq r))))
    (assoc m :blast2-db db)))

(defn- results-db
  [{:keys [wdir] :as m}]
  (let [db (bdb/db-spec {:dbname (str (fs/file wdir "rec-blast-results.sqlite"))
                         :dbtype :sqlite})]
    (bdb/create-table! db :sequences :blast)
    (assoc m :results db)))

(defn filter-save-reciprocal-hits
  [{:keys [blast1 blast2-db results relax] :as m}]
  (mapcat #(with-open [r (io/reader %)]
             (let [res (->> (bl/iteration-seq r)
                            (filter (fn [x] (blasts-match? x blast2-db relax))))]
               (bdb/insert-sequences! results :sequences :blast res)))
          blast1)
  m)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; api
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn reciprocal-blast
  "Runs reciprocal blast on a collection of fasta sequences."
  ([coll database program] (reciprocal-blast coll database program {}))
  ([coll database program {:keys [wd evalue bparams relax query-coverage]
                           :or {wd nil evalue 10e-5 bparams {} relax 1 query-coverage 0.70}}]
   (let [wdir (make-wd wd)
         m {:wdir wdir
            :sfile (fa/fasta->file coll (fs/file wdir "sequences1.txt"))
            :database database
            :relax relax
            :program program
            :query-coverage query-coverage
            :evalue evalue
            :bparams (merge bparams {"-max_target_seqs" (str relax)})}]
     (->> (first-blasts m)
          second-blasts
          save-blasts
          results-db
          filter-save-reciprocal-hits))))

(defn reciprocal-blast-file
  "Runs reciprocal blast on a file of fasta formatted sequences."
  ([file database program] (reciprocal-blast-file file database program {}))
  ([file database program {:keys [wd evalue bparams relax query-coverage]
                           :or {wd nil evalue 10e-5 bparams {} relax 1 query-coverage 0.70}}]
   (let [wdir (make-wd wd)
         m {:wdir wdir
            :sfile file
            :relax relax
            :database database
            :program program
            :query-coverage query-coverage
            :evalue evalue
            :bparams (merge bparams {"-max_target_seqs" (str relax)})}]
     (->> (first-blasts m)
          second-blasts
          save-blasts
          results-db
          filter-save-reciprocal-hits))))

