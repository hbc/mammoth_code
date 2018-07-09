import sys
import os

import subprocess


out = "%s-snpeff.tsv" % (os.path.splitext(sys.argv[1])[0])

with open(sys.argv[1]) as inh:
    with open(out, 'w') as outh:
        counts = 0
        for line in inh:
            counts += 1
            if counts % 1000 == 0:
                print "reading %s hits" % counts
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.split()[9:]
                print >>outh, "chrom\tpos\tref\talt\ttype\timpact\tgene\tprotein\twarning"
                continue
            cols = line.split("\t")
            chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
            ann = cols[7]
            snpeff = [atr for atr in ann.split(";") if atr.startswith("ANN")]
            if not snpeff[0]:
                continue
            snpeff = snpeff[0].split("|")
            if snpeff[1].find("missense") > -1 or snpeff[1].find("stop") > -1:
                print >>outh, "\t".join([chrom, pos, ref, alt, snpeff[1], snpeff[2], snpeff[3], snpeff[10], snpeff[15]])
