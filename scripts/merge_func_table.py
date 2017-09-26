import os
import sys
import pybedtools
import shutil
from collections import defaultdict

hg19="ann/hg19_pos.bed"
variants="all_23genomes.xls"
ann="ann/hg19_ann.txt"

hg19_map = {}
with open(hg19) as inh:
    for line in inh:
        cols = line.strip().split()
        hg19_map[cols[3]] = "%s:%s" % (cols[0], int(cols[1]) + 1)

ann_map = defaultdict(dict)
with open(ann) as inh:
    for line in inh:
        cols = line.strip().split()
        idx = "chr%s:%s" % (cols[7], cols[8])
        cols[4] = cols[4] if cols[4] != "X" else "STOP"
        cols[5] = cols[5] if cols[5] != "X" else "STOP"
        ann_map[idx].update({"%s %s" % (cols[4],cols[5]): cols[34]})
        # print ann_map[idx]

def _get_change(cols):
    if cols[7] != "NA":
        ref = cols[7]
    elif cols[27]:
        ref = aa_name[cols[27][2:5]]
    if cols[11] != "NA":
        alt = cols[11]
    elif cols[18] != "NA":
        alt = cols[18]
    else:
        if cols[27].endswith("*"):
            alt = "STOP"
        elif cols[27][-3:] not in aa_name:
            alt = "NA"
        else:
            alt = aa_name[cols[27][-3:]]
    return [ref, alt]

aa_name = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
           'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
           'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
           'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
           'Ter': 'STOP', '': 'NA'}

with open(variants) as inh:
    for line in inh:
        if line.find("gene") > -1:
            print "\t".join(line.strip().split("\t") + ["polyphen2"])
            continue
        cols=[c.replace("\"", "") for c in line.strip().split("\t")]
        chrom, pos = cols[0], int(cols[1])
        idx = "%s:%s" % (chrom.replace("\"", ""), pos)
        score = "NA"
        if idx in hg19_map:
            if hg19_map[idx] in ann_map:
                change = _get_change(cols)
                score = []
                for item in ann_map[hg19_map[idx]]:
                    score.append(ann_map[hg19_map[idx]][item])
                    if " ".join(change) == item:
                        score = score[-1]
                        break
                score = ",".join(score)
        print "\t".join(cols + [score])
