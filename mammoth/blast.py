import json
import yaml
import os
from collections import defaultdict

import mammoth.logger as mylog
from mammoth import utils, do
from mammoth.ensembl import query_exon, query_prot

logger = mylog.getLogger()

def _read_json(st):
    jstxt = ""
    res = []
    for line in st.split("\n"):
        if line == "}":
            jstxt += line
            res.append(json.loads(jstxt))
            jstxt = ""
        else:
            jstxt += line
    return res

def find(s, ch = "X"):
    """http://stackoverflow.com/a/11122355/1772223"""
    return [i for i, ltr in enumerate(s) if ltr == ch]

def _identify_mut(qseq, matches, hit, start, tx, prot):
    mism = find(matches)
    res = dict()
    if len(mism) > 10:
        return {'positions': res}
    for m in mism:
        abs_pos = start + m
        aa = int(abs_pos / 3) + 1
        logger.debug("mism found %s at pos %s" % (m, abs_pos))
        logger.debug("aa at: %s" % aa)
        shift = (aa - 1) * 3 - abs_pos
        logger.debug("shift %s" % shift)
        logger.debug("seq %s" % qseq)
        logger.debug("prot %s aas" % len(prot['seq']))
        if aa > len(prot['seq']):
            continue
        new_aa = triplet.get(qseq[(m+shift):(m+shift+3)], None)
        ref_aa = prot['seq'][aa-1]
        ref_tri = tx['seq'][(aa-1)*3:aa*3]
        new_tri = qseq[(m+shift):(m+shift+3)]
        logger.debug("triplete %s->%s with aa %s->%s at %s" % (ref_tri, new_tri, ref_aa, new_aa, aa))
        if ref_aa != new_aa:
            res[m]={'ref_nc': ref_tri, 'ref_aa': ref_aa, 'new_nc': new_tri, 'new_aa': new_aa, 'aa_pos': aa, 'tx_pos': abs_pos, 'mism_pos': m}
    return {'positions': res}

def _consistency(res, gene):
    exon_sort = dict(sorted(gene['exons'].iteritems()))
    pos = []
    hits = []
    qc = defaultdict(dict)
    chrom = None
    for e in exon_sort:
        logger.debug(exon_sort[e])
        name = exon_sort[e]
        if name in res:
            if not chrom:
                chrom = res[name]['chr']
            res[name]['exon_number'] = e
            if chrom == res[name]['chr']:
                pos.append(res[name]['pos'])
                hits.append(e)
            else:
                qc['other_chr'].update({e: name})
        else:
            qc['missing'].update({e: name})
    qc.update({'mapped': len(pos), 'annotated': len(gene['exons'])})

    if len(pos) > 0:
        size = abs(min(pos) - max(pos))
        logger.debug(pos)
        logger.debug([size, gene['size'] * 2])
        logger.debug(gene['size'] * 2 > size)
        if gene['size'] * 2 > size:
            return res, qc
    return None, qc

def _parse_algn(algn, tx, prot, gene):
    res = dict()
    for a in algn:
        r = a['BlastOutput2']
        if 'report' in r:
            hits =  r['report']['results']['search']['hits']
            name = r['report']['results']['search']['query_title']
            for h in hits:
                chrom = h['description'][0]['id']
                if name in res:
                    break
                if chrom.find("Un") < 0:
                    qseq = str(h['hsps'][0]['qseq']).upper()
                    hseq = str(h['hsps'][0]['hseq']).upper()
                    pos_hit = h['hsps'][0]['hit_from']
                    matches = str(h['hsps'][0]['midline']).replace(" ", "X")
                    pos_start = tx["seq"].find(qseq)
                    logger.debug(h)
                    res.update({name: {'chr': chrom,'pos': pos_hit, 'matches' : matches,
                                       'pos_start': pos_start, 'hseq': hseq, 'qseq': qseq}})
    res, qc = _consistency(res, gene)
    if res:
        logger.debug("Matches %s out of %s" % (qc['mapped'], qc['annotated']))
        gap = 0
        for name in res:
            logger.debug(res[name])
            if res[name]['qseq'].find("-") < 0 and res[name]['hseq'].find("-") < 0:
                res[name]['changes'] = _identify_mut(res[name]['hseq'], res[name]['matches'],
                                                     res[name]['pos'], res[name]['pos_start'],
                                                     tx, prot)
            else:
                gap = 1
        qc.update({'gap': gap})
        res['mapped_exons'] = qc
    return res

def blast_genes(genes, db, outdir):
    outdir = utils.safe_makedir(outdir)
    cmd = "blastn -query {fn} -db {db} -task blastn -qcov_hsp_perc 90 -outfmt 13 > {out}"
    for g in genes:
        for t in genes[g]:
            fn = os.path.join(outdir, "%s.input" % t)
            out = "%s.blastn" % fn
            out_final = "%s.changes" % fn
            if utils.file_exists(out_final):
                changes = yaml.load(open(out_final))
                genes[g][t]['changes'] = changes
                continue
            tx_seq = query_exon(t)
            protein = query_prot(t)
            exon_sort = dict(sorted(genes[g][t]['exons'].iteritems()))
            if not utils.file_exists(fn):
                with open(fn, 'w') as ih:
                    for e in exon_sort:
                        seq = query_exon(exon_sort[e])
                        logger.debug(seq)
                        print >>ih, ">%s\n%s" % (seq["id"], seq["seq"])
            if not utils.file_exists(out):
                do.run(cmd.format(**locals()))
            if not utils.file_exists(out_final):
                algn = _read_json(open(out).read())
                changes = _parse_algn(algn, tx_seq, protein, genes[g][t])
                with open(out_final, 'w') as outh:
                    yaml.dump(changes, outh, default_flow_style=False, allow_unicode=False)
            genes[g][t]['changes'] = changes
    return genes

triplet = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
"TCT":"S","TCC":"s","TCA":"S","TCG":"S",
"TAT":"Y","TAC":"Y","TAA":"STOP","TAG":"STOP",
"TGT":"C","TGC":"C","TGA":"STOP","TGG":"W",
"CTT":"L","CTC":"L","CTA":"L","CTG":"L",
"CCT":"P","CCC":"P","CCA":"P","CCG":"P",
"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
"CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M",
"ACT":"T","ACC":"T","ACA":"T","ACG":"T",
"AAT":"N","AAC":"N","AAA":"K","AAG":"K",
"AGT":"S","AGC":"S","AGA":"R","AGG":"R",
"GTT":"V","GTC":"V","GTA":"V","GTG":"V",
"GCT":"A","GCC":"A","GCA":"A","GCG":"A",
"GAT":"D","GAC":"D","GAA":"E","GAG":"E",
"GGT":"G","GGC":"G","GGA":"G","GGG":"G",}
