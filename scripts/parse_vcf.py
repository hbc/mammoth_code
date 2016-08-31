import sys
import os

import subprocess

def get_genotype(samples):
    gen = []
    for s in samples:
        var = s.split(":")[0]
        if var == "1/0":
            gen.append("Het")
        elif var == "1/1":
            gen.append("Hom")
        else:
            gen.append("None")
    return gen

def get_seq(chrom, pos, bam):
    """use samtools to get consensus"""
    tmp = "tmp.samtools"
    cmd = "samtools tview -d T -p {0}:{1} {2}".format(chrom, pos, bam)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    # process.communicate()
    inh = process.stdout
    # process.stdout.close()
    pos = inh.readline()
    ref = inh.readline()
    seq = inh.readline().strip()
    inh.close()
    return seq

def get_flank(chrom, pos, names, gen):
    flank_list = []
    for i, sample in enumerate(names):
        if gen[i] != "Hom":
            flank_list.append("None")
            continue
        bam = os.path.join("~/scratch/church_mammoth/mammoth_vc/work/align", sample, "%s-sort.bam" % sample)
        pre = int(pos) - 80
        first = get_seq(chrom, pre, bam)
        second = get_seq(chrom, pos, bam)
        flank = "NoInfo"
        if len(second) > 2:
            flank = "%s-%s-%s" % (first.replace(" ", "X"), second[0], second[1:].replace(" ", "X"))
        flank_list.append(flank)
    return flank_list

def _get_cache(fn):
    cache = {}
    with open(fn) as inh:
        for line in inh:
            line = line.strip()
            cols = line.split()
            idx = "%s%s" % (cols[0], cols[1])
            cache[idx] = line
    return cache

out = "%s-parsed-flank.tsv" % (os.path.splitext(sys.argv[1])[0])

cache = {}
if os.path.exists(out):
    cache = _get_cache(out)

with open(sys.argv[1]) as inh:
    with open(out, 'w') as outh:
        counts = 0
        for line in inh:
            counts += 1
            if counts % 100 == 0:
                print "reading %s hits" % counts
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.split()[9:]
                print >>outh, "chrom\tpos\tref\talt\tgene\tchange\t%s\t%s" % ("\t".join(samples), "\t".join(samples))
                continue
            cols = line.split("\t")
            chrom, pos, alt, ref = cols[0], cols[1], cols[3], cols[4]
            ann = cols[7]
            snpeff = ann.split(";")[-1].split("|")
            if snpeff[1].find("missense") > -1:
                idx = "%s%s" % (chrom, pos)
                if idx in cache:
                    print >>outh, cache[idx]
                else:
                    # print [chrom, pos]
                    gen = get_genotype(cols[9:])
                    flanks = get_flank(chrom, pos, samples, gen)
                    print >>outh, "\t".join([chrom, pos, alt, ref, snpeff[3], snpeff[10]] + gen + flanks)
