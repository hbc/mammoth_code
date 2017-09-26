export PATH=~/scratch/church_mammoth/conda/bin:$PATH

mammoth-run.py annotate --gtf Loxodonta_africana.loxAfr3.85.gtf --db ~/mammoth/Oimyakon/OimyakonMammothBlastDb --fasta ~/mammoth/Oimyakon/Oimyakon_consensus_sequence.fa -o res_oimyako -n 2 -d --fasta_ref /groups/bcbio/bcbio/genomes/Lafricana/loxAfr3/seq/loxAfr3.fa
