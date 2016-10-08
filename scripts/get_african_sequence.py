import os
import sys
import pybedtools
import shutil

genome="/groups/bcbio/bcbio/genomes/Lafricana/loxAfr3/seq/loxAfr3.fa"
size="/groups/bcbio/bcbio/genomes/Lafricana/loxAfr3/seq/loxAfr3.fa.fai"
fai = dict()
with open(size) as inh:
    for line in inh:
        cols=line.strip().split()
        chrom, size = cols[0], int(cols[1])
        fai[chrom] = size

variants="/home/lp113/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1/split/merged-parsed-flank-wheader.tsv"
# variants="../test/data/test.vcf"
def _get_cache(fn):
    cache = {}
    with open(fn) as inh:
        for line in inh:
            line = line.strip()
            cols = line.split()
            idx = "%s%s" % (cols[0], cols[1])
            cache[idx] = line
    return cache

out_file = "merged.fa"
tmp_file = "merged.tmp.fa"
cache = {}
if os.path.exists(out_file):
    cache = _get_cache(out_file)
    shutil.move(out_file, tmp_file)

with open(out_file, "w") as outh:
    with open(variants) as inh:
        for line in inh:
            if line.startswith("chrom"):
                continue
            cols=line.strip().split()
            chrom, pos=cols[0], int(cols[1])
            idx = "%s%s" % (chrom, pos)
            if idx in cache:
                print >>outh, cache[idx]
                continue

            start = 0 if pos - 150 < 0 else pos - 150
            end = fai[chrom] - 1 if pos + 150 - 1 > fai[chrom] else pos + 150 - 1
            a = pybedtools.BedTool("{0}\t{1}\t{2}".format(chrom, start, end), from_string=True)
            fasta = pybedtools.example_filename(genome)
            a = a.sequence(fi=fasta)
            seq = open(a.seqfn).read().split("\n")
            pre = seq[1][:pos-start-1]
            nt = seq[1][pos-start-1]
            post = seq[1][pos-start:]
            # print [pre, nt, pos]
            print >>outh, "%s\t%s\t%s-%s-%s" % (chrom, pos, pre, nt, post)
