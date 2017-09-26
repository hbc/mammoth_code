export PATHROOT=/home/lp113/temporal/lp113/mammoth/elephants
export VCF=b1-gatk-haplotype-joint-annotated-decomposed.vcf

cd elephants/final/2017-09-06_elephants
mkdir split && cd split
mkdir log

gunzip ../$VCF.gz

split ../$VCF

grep "^#" ../$VCF > header

for F in `ls x* | grep -v tsv | grep -v "\.h"` ; do
    echo $F
    `cat header $F > $F.h`
    if [ ! -e $F"-parsed-flank.tsv"  ] ; then
        # echo do
        /usr/bin/srun -o log/$F.o -e log/$F.e -t 0-13:00 --mem 8000 -n 1 -p medium  $PATHROOT/conda/bin/python $PATHROOT/mammoth_code/scripts/parse_vcf.py $F.h &
    fi
done

rm -rf  merged-parsed-flank.tsv
grep "^scaffold" *-parsed-flank.tsv | sed 's/^.*-parsed-flank.tsv://'  > merged-parsed-flank.tsv

cat <(grep -w chrom xab-parsed-flank.tsv)  <(cat merged-parsed-flank.tsv) > merged-parsed-flank-wheader.tsv

awk '{start=$2-150; if (start<0){start=0}; print $1"\t"$2-150"\t"$2+150}' ../$VCF |grep -v "##" > merged.bed

/usr/bin/srun -o log/fasta.o -e log/fasta.e --mem 8000 -t 0-13:00 -n  1 -p medium  $PATHROOT/conda/bin/python $PATHROOT/mammoth_code/scripts/get_african_sequence.py  merged-parsed-flank-wheader.tsv

# convert to human to use snpsniff
grep -v "^#"  ../$VCF | awk '{print $1"\t"$2"\t"$2+1"\t"$1":"$2}' > vcf_pos.bed

# get exact genotype
/usr/bin/srun -o log/genotype.o -e log/genotype.e -t 0-40:00 --mem 8000 -n  1 -p medium  $PATHROOT/conda/bin/python $PATHROOT/mammoth_code/scripts/parse_vcf_genotype.py ../$VCF

cd $PATHROOT/mammoth_code/table
Rscript ../scripts/merge-tables2.R
bash all_genome_ann.sh
Rscript ../scripts/clean-tables.R

