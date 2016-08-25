from joblib import Parallel, delayed
import sys
import os
import time

from mammoth.logger import initialize_logger
from mammoth.parse import parse_cl
import mammoth.logger as mylog
from mammoth import ensembl
from mammoth.blast import blast_genes

def _get_genomic(m, start, end, strand):
    if strand=="-":
        return end - m
    else:
        return start + m

def _get_genomic_sequence(exons, pos):
    seq = ensembl.query_sequence(exons['chrom'], pos-150, pos+150, exons['strand'])
    return "%s-%s-%s" % (seq[:150], seq[150], seq[151:])

def _format(info, exons):
    genomic_pos = _get_genomic(info[3]['mism_pos'], exons['start'], exons['end'], exons['strand'])
    flank = _get_genomic_sequence(exons, genomic_pos)
    changes = "%d %s %s %d %s %s" % (info[3]['tx_pos'], info[3]['ref_nc'], info[3]['new_nc'],
                                     info[3]['aa_pos'], info[3]['ref_aa'], info[3]['new_aa'])
    missing = ",".join(["%s:%s" %  (e, info[4]['missing'][e]) for e in info[4]['missing']])
    qc = "%s/%s %s %s" % (info[4]['mapped'], info[4]['annotated'], info[4]['gap'], missing)
    return "%s %s %s %s %s %s %s %s %s" % (info[0], info[1], info[2], info[5], changes, exons['chrom'], genomic_pos, qc, flank)

def _chunk(l, n):
    i = 0
    idx = []
    jump = int(1.0 * l / n)
    while len(idx) < n:
        idx.append([i, i + jump])
        i += jump
    idx[-1][1] = l
    return idx

def _join(matches):
    out = dict()
    for i in matches:
        for g in i:
            out.update({g: i[g]})
    return out

def _get_cache(fn):
    cache = {}
    with open(fn) as inh:
        for line in inh:
            line = line.strip()
            cols = line.split()
            idx = "%s%s%s%s" % (cols[0], cols[1], cols[2], cols[4])
            cache[idx] = line
    return cache

def run_smartly(genes, args):
    logger.info("Runnning %s genes" % len(genes.keys()))
    if args.n == 1:
        matches = blast_genes(genes, args.db, args.out)
        return matches
    idx = _chunk(len(genes.keys()) , args.n)
    list_genes = []
    last = 0
    for i in idx:
        out = dict()
        [out.update({k: genes[k]}) for k in genes.keys()[i[0]:i[1]]]
        list_genes.append(out)
    logger.info("Running %s genes in %s chunks" % (len(genes.keys()), len(list_genes)))
    matches = Parallel(args.n)(delayed(blast_genes)(out, args.db, args.out) for out in list_genes)
    return _join(matches)

def run(args):
    """Proxy for the pipeline"""
    genes, exons = ensembl.get_genes(args.gtf)
    matches = run_smartly(genes, args)
    out_file = os.path.join(args.out, "changes.tsv")
    if os.path.exists(out_file):
        cache = _get_cache(out_file)
    with open(out_file, 'w') as outh:
        print >>outh, " ".join(["gene", "tx", "exon", "number", "tx_pos", "ref_nc", "ref_new",
                                "aa_pos", "ref_aa", "new_aa", "chrom", "genome_pos", "mapped/annotated",
                                "found_gap", "missing_exons", "flank"])
        for m in matches:
            for t in matches[m]:
                if not matches[m][t]['changes']:
                    continue
                ratio = matches[m][t]['changes']['mapped_exons']
                for e in matches[m][t]['changes']:
                    if "changes" in matches[m][t]['changes'][e]:
                        for p in matches[m][t]['changes'][e]['changes']['positions']:
                            idx = "%s%s%s%s" % (m, t, e, p)
                            if idx in cache:
                                print >>outh, cache[idx]
                                continue
                            info = matches[m][t]['changes'][e]['changes']['positions'][p]
                            en = matches[m][t]['changes'][e]['exon_number']
                            print >>outh, _format([m, t, e, info, ratio, en], exons[e])

if __name__ == "__main__":
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out, kwargs['args'].debug, kwargs['args'].print_debug)
    logger = mylog.getLogger()
    start = time.time()
    if "annotate" in kwargs:
        logger.info("Run annotation")
        run(kwargs["args"])
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))


