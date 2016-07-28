import argparse
import sys


def parse_cl(in_args):
    print in_args
    sub_cmds = {"annotate": add_subparser_mirbuster,
                }
    parser = argparse.ArgumentParser(description="small RNA analysis")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="seqcluster supplemental commands")
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    else:
        print "use %s" % sub_cmds.keys()
        sys.exit(0)
    args = parser.parse_args()

    assert sub_cmd is not None
    kwargs = {"args": args, sub_cmd: True}
    return kwargs


def _add_debug_option(parser):
    parser.add_argument("-d", "--debug", action="store_true",
                        dest="debug", help="max verbosity mode", default=False)
    parser.add_argument("-vd", "--print_debug", action="store_true",
                        help="print debug messageson terminal", default=False)
    return parser


def add_subparser_mirbuster(subparsers):
    parser = subparsers.add_parser("annotate", help="realign miRNA BAM file")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--gtf", help="gtf file with precursor position to genome.", required=1)
    parser.add_argument("--db", help="Blast database.", required=1)
    parser.add_argument("-n", help="Jobs.", type=int, default=1)
    parser = _add_debug_option(parser)
    return parser


