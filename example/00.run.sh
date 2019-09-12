perl ../bin/wgddetector.pl --input_cds test.cds.fa --input_pep test.pep.fa --output_dir out --tmp_dir tmp --thread_num 4 --cluster_engine mmseqs2 --clean no

Rscript example.Rscript

ln -s output/04.final_paralogs_ks/final.ks.distribution.list.pdf .
