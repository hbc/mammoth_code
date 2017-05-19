# argument to be fiving to bash script: ~/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1
cd $1
SOURCE=`dirname $0`
PYTHON="~/scratch/church_mammoth/conda/bin/"

mkdir split && cd split
mkdir log

gunzip ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf.gz

split ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf

grep "^#" ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf > header

for F in `ls x* | grep -v tsv | grep -v "\.h"` ; do
    echo $F
    `cat header $F > $F.h`
    if [ ! -e $F"-parsed-flank.tsv"  ] ; then
        # echo do
        bsub -oo log/$F.o -W 10:00 -R "rusage[mem=8000]" -n  2 -q parallel  $PYTHON/python $SOURCE/parse_vcf.py $F.h
    fi
done

rm -rf  merged-parsed-flank.tsv
grep "^scaffold" *-parsed-flank.tsv | sed 's/^.*-parsed-flank.tsv://'  > merged-parsed-flank.tsv

cat <(grep -w chrom xab-parsed-flank.tsv)  <( cat merged-parsed-flank.tsv) > merged-parsed-flank-wheader.tsv

awk '{start=$2-150; if (start<0){start=0}; print $1"\t"$2-150"\t"$2+150}' ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf |grep -v "##" > merged.bed

bsub -eo log/fasta.o -W 40:00 -R "rusage[mem=8000]" -n  2 -q parallel  $PYTHON/python $SOURCE/get_african_sequence.py

# convert to human to use snpsniff
grep -v "^#"  ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf | awk '{print $1"\t"$2"\t"$2+1"\t"$1":"$2}' > vcf_pos.bed

# get exact genotype
bsub -eo log/genotype.o -W 40:00 -R "rusage[mem=8000]" -n  1 -q priority  $PYTHON/python $SOURCE/parse_vcf_genotype.py ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf

cd $SOURCE/../table
Rscript ../scripts/merge-tables.R
bash all_genome_ann.sh
Rscript ../scripts/clean-tables.R

