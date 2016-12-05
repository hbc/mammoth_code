import sys
import os

import subprocess

def get_genotype(samples, ref, alt):
    gen = []
    nt = [ref] + alt.split(",")
    for s in samples:
        var = s.split(":")[0].split("/")
        g1 = '.' if var[0] == '.' else nt[int(var[0])]
        g2 = ''
        if len(var) == 2:
            g2 = '.' if var[1] == '.' else nt[int(var[1])]
        gen.append("%s/%s" % (g1,g2))
    return gen

def _get_cache(fn):
    cache = {}
    with open(fn) as inh:
        for line in inh:
            line = line.strip()
            cols = line.split()
            idx = "%s%s" % (cols[0], cols[1])
            cache[idx] = line
    return cache

out = "%s-parsed-genotype.tsv" % (os.path.splitext(sys.argv[1])[0])

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
                print >>outh, "chrom\tpos\t%s" % ("\t".join(samples))
                continue
            cols = line.split("\t")
            chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
            ann = cols[7]
            snpeff = [atr for atr in ann.split(";") if atr.startswith("ANN")]
            if not snpeff[0]:
                continue
            snpeff = snpeff[0].split("|")
            if snpeff[1].find("missense") > -1 or snpeff[1].find("stop") > -1:
                idx = "%s%s" % (chrom, pos)
                # print [chrom, pos]
                try:
                    gen = get_genotype(cols[9:], ref, alt)
                    print >>outh, "\t".join([chrom, pos] + gen)
                except:
                    print line
                    raise IndexError("error in genotype format")
