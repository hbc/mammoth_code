from unittest import TestCase
import os
import inspect
from nose.plugins.attrib import attr

import logging
from mammoth import blast
from mammoth.logger import initialize_logger
import mammoth.logger as mylog

initialize_logger(None, True, True)
logger = mylog.getLogger(__name__)

class TestFunctions(TestCase):
    @attr(mutation=True)
    def test_alignment(self):
        tx = {"id":"ENSLAFT00000013206'","seq":"ATGAATGCTCACCCCAAGGAGATGGTGCCCCTTATGGGCAGGAAAGCCATGGCCCCCAGTGGGAACCCTGCCGTCCTACAGGAGAAGAGGCCTGCAGAGCTCACCCCTACCAAGAAAAG"}
        prot={"id":"ENSLAFT00000013206'","seq":"MNAHPKEMVPLMGRKAMAPSGNPAVLQEKRPAELTPTKK"}
        gene={"exons":{1:"exon1"}, "size":1000}
        out = os.path.join(os.path.dirname(__file__), "align.blast")
        a = blast._read_json(open(out).read())
        blast._parse_algn(a, tx, prot, gene)
