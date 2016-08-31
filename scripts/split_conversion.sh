
for F in `ls * | grep -v tsv | grep -v "h$"` ; do
    echo $F
    `cat header $F > $F.h`
    if [ ! -e $F"-parsed-flank.tsv"  ] ; then
        # echo do
        bsub -W 10:00 -R "rusage[mem=8000]" -n  2 -q mcore  ~/scratch/church_mammoth/conda/bin/python ~/scratch/church_mammoth/mammoth_code/scripts/parse_vcf.py $F.h
    fi
done

rm -rf  merged-parsed-flank.tsv
grep "^scaffold" *-parsed-flank.tsv | sed 's/^.*-parsed-flank.tsv://'  > merged-parsed-flank.tsv

cat <(grep -w chrom xab-parsed-flank.tsv)  <( cat merged-parsed-flank.tsv) > merged-parsed-flank-wheader.tsv

awk '{start=$2-150; if (start<0){start=0}; print $1"\t"$2-150"\t"$2+150}' ../batch1-joint-effects-filterSNP-filterINDEL-gatkclean.vcf |grep -v "##" > merged.bed

betdo
