import sys

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

with open(sys.argv[1]) as inh:
    for line in inh:
        line = line.strip()
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.split()[9:]
            print "chrom\tpos\tref\talt\tgene\tchange\t%s" % "\t".join(samples)
            continue
        cols = line.split("\t")
        chrom, pos, alt, ref = cols[0], cols[1], cols[3], cols[4]
        ann = cols[7]
        snpeff = ann.split(";")[-1].split("|")
        if snpeff[1].find("missense") > -1:
            gen = get_genotype(cols[9:])
            print "\t".join([chrom, pos, alt, ref, snpeff[3], snpeff[10]] + gen)
