"""ensembl interaction function"""
import os
import requests, sys
import yaml
import logging

import gffutils

from collections import defaultdict

import mammoth.logger as mylog

server = "http://rest.ensembl.org{ext}"
ext = "/sequence/id/{id}?type=cds"
prot = "/sequence/id/{id}?type=protein"
sequence = "/sequence/region/elephant/{chr}:{start}..{end}:{strand}?"

def query_sequence(chr, start, end, strand):
    r = requests.get(server.format(ext=sequence.format(**locals())), headers={ "Content-Type" : "text/plain"})
    if not r.ok:
        r.raise_for_status()
        return None
    return yaml.load(r.text)

def query_exon(id):
    r = requests.get(server.format(ext=ext.format(id=id)), headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        return None
    return yaml.load(r.text)

def query_prot(id):
    r = requests.get(server.format(ext=prot.format(id=id)), headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        return None
    return yaml.load(r.text)

def _get_db(db):
    return gffutils.FeatureDB(db_file)

def _convert_to_db(db):
    out = "%s.db" % db
    if os.path.exists(out):
        return gffutils.FeatureDB(out)
    gffutils.create_db(db, disable_infer_transcripts=True, disable_infer_genes=True, dbfn=out)
    return gffutils.FeatureDB(out)

def get_genes(db):
    db = _convert_to_db(db)
    genome = defaultdict(dict)
    exons_pos = defaultdict(dict)
    for gene in db.features_of_type("gene"):
        if "gene_name" not in gene.attributes:
            continue
        if gene.attributes["gene_biotype"][0] == "protein_coding":
            exon_seen = set()
            for tx in db.children(gene, featuretype='transcript', order_by='start'):
                if tx.attributes["transcript_biotype"][0] == "protein_coding":
                    # txs.add(tx["transcript_id"])
                    exons = dict()
                    for e in db.children(tx, featuretype='exon', order_by='start'):
                        if e.attributes['exon_id'][0] not in exon_seen:
                            exons.update({int(e.attributes['exon_number'][0]): e.attributes['exon_id'][0]})
                            exons_pos.update({e.attributes['exon_id'][0]: {'chrom': e.chrom,
                                                                           'start': e.start,
                                                                           'end': e.end,
                                                                           'strand': e.strand}})
                        exon_seen.add(e.attributes['exon_id'][0])
                    genome[gene.attributes["gene_name"][0]].update({tx.attributes["transcript_id"][0]: {'size': abs(tx.end-tx.start),
                    'exons': exons}})
    return genome, exons_pos

