# clj-reciprocal-blast

Clojure library for running reciprocal blast.

## Installation

To include in project: [clj-reciprocal-blast "0.1.2"]

To include namespace: (:require [clj-reciprocal-blast.core :as rb])

## Usage

Provides two functions `reciprocal-blast` and `reciprocal-blast-file`
that take a blast database, a program and either a collection of fasta
sequences (see `clj-fasta`) or a file of fasta formatted sequences,
respectively, and runs reciprocal blast on the sequences. Additional
parameters (working directory, evalue cutoff, blast parameters and the
number of reblasted hits to take into account for relaxed blast) can
be supplied as an optional hashmap. Returns a hashmap containing an
sqlite database with the results as well as the location of all the
intermediate files:

```clojure
user> (with-open [r (io/reader "/path/to/fasta-seqs.fasta")]
        (reciprocal-blast (fa/fasta-seq r) "/path/blast-db" "blastp"))
{:wdir "./wdir/", :sign-fasta "./wdir/bl1-significant-hits.fasta",
 :relax 1, :sfile #object[java.io.File 0x5e68ea05 "./wdir/sequences1.txt"],
 :blast2 ("./wdir/blast2-1.xml-1.xml"), :bparams {"-max_target_seqs" "1"},
 :blast1 ("./wdir/blast1-1.xml-1.xml"), :database "/path/blast-db",
 :program "blastp", :blast2-db {:classname "org.sqlite.JDBC",
 :subprotocol "sqlite", :subname "./wdir/blast2.sqlite", :dbtype :sqlite},
 :evalue 1.0E-4, :results {:classname "org.sqlite.JDBC", :subprotocol "sqlite",
 :subname "./wdir/rec-blast-results.sqlite", :dbtype :sqlite}}
user>
```

The value of :results in the returned hashmap is a jdbc database spec
refering to a sqlite file produced at the end of the analysis that
contains blast results for each original sequence that produced a
reciprocal hit. Results can be accessed using `biodb` and `clj-blast`:

```clojure
user> (bdb/query-sequences (:results rblast-res) ["select * from sequences"] :blast
                           :apply-func #(doseq [r (->> (map bl/hit-seq %)
                                                       (map first))]
                                          (println (:Hit_def x))))
Neoverrucotoxin subunit beta OS=Synanceia verrucosa PE=1 SV=1
Calglandulin OS=Bothrops insularis PE=2 SV=1
Cysteine-rich venom protein 1 OS=Pimpla hypochondriaca PE=1 SV=1
C-type lectin mannose-binding isoform OS=Notechis scutatus scutatus PE=2 SV=1
U3-aranetoxin-Ce1a OS=Caerostris extrusa PE=2 SV=1
...
user>
```

## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
