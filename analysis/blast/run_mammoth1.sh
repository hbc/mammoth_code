export PATH=~/scratch/church_mammoth/conda/bin:$PATH

mammoth-run.py annotate --gtf Loxodonta_africana.loxAfr3.85.gtf --db ~/mammoth/Wrangel/WrangelMammothBlastDb  --fasta ~/mammoth/Wrangel/Wrangel_consensus_sequence.fa  -o res_wrangel -n 2 -d --fasta_ref /groups/bcbio/bcbio/genomes/Lafricana/loxAfr3/seq/loxAfr3.fa
