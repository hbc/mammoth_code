cd ~/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1
mkdir split && cd split
mkdir log

gunzip ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf.gz

split ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf

grep "^#" ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf > header

for F in `ls x* | grep -v tsv | grep -v ".h"` ; do
    echo $F
    `cat header $F > $F.h`
    if [ ! -e $F"-parsed-flank.tsv"  ] ; then
        # echo do
        bsub -oo log/$F.o -W 10:00 -R "rusage[mem=8000]" -n  2 -q parallel  ~/scratch/church_mammoth/conda/bin/python ~/scratch/church_mammoth/mammoth_code/scripts/parse_vcf.py $F.h
    fi
done

rm -rf  merged-parsed-flank.tsv
grep "^scaffold" *-parsed-flank.tsv | sed 's/^.*-parsed-flank.tsv://'  > merged-parsed-flank.tsv

cat <(grep -w chrom xab-parsed-flank.tsv)  <( cat merged-parsed-flank.tsv) > merged-parsed-flank-wheader.tsv

awk '{start=$2-150; if (start<0){start=0}; print $1"\t"$2-150"\t"$2+150}' ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf |grep -v "##" > merged.bed

bsub -eo log/fasta.o -W 40:00 -R "rusage[mem=8000]" -n  2 -q parallel  ~/scratch/church_mammoth/conda/bin/python ~/scratch/church_mammoth/mammoth_code/scripts/get_african_sequence.py

# conver to human to use snpsniff
grep -v "^#"  ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf | awk '{print $1"\t"$2"\t"$2+1"\t"$1":"$2}' > vcf_pos.bed

# CrossMap.py bed ~/scratch/church_mammoth/tools/liftover/loxAfr3ToHg38.over.chain.gz vcf_pos.bed vcf_pos_hg38.bed

# ~/scratch/church_mammoth/conda/bin/python ~/scratch/church_mammoth/mammoth_code/scripts/merge_bed_vcf.py vcf_pos_hg38.vcf
# bedtools sort -i vcf_pos_hg38.vcf.gz -faidx /groups/bcbio/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa.fai -header > vcf_pos_sort_hg38.vcf
# bgzip vcf_pos_hg38_sort.vcf
# bcftools index vcf_pos_hg38_sort.vcf.gz
# java -jar ~/scratch/church_mammoth/tools/snpEff/SnpSift.jar dbnsfp -db ~/scratch/church_mammoth/tools/snpEff/db/hg19/dbNSFP/dbNSFP3.2.txt.gz  vcf_pos_hg38_sort.vcf.gz


