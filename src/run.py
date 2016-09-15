import os
import sys
import argparse
import time
import re
import utils
import constants
from collections import defaultdict as dd
from jinja2 import Environment, FileSystemLoader, escape
import json

EXEC_DIR = utils.module_path() + "/"


def read_exons(exonsf):
    res = {}
    for l in exonsf:
        chrom, strand, gene, exons = l.rstrip().split('\t')
        if '(' in gene:  # TODO: temporary solution for current files -- REMOVE
            gene = re.match(".*\((.*)\).*", gene).groups()[0]
        res[gene] = {'chrom': chrom, 'strand': strand, 'exons': exons}
    return res


def read_mat(mat):
    fields = {}
    for l in mat:
        if l.startswith('gene_id'): continue
        k, v = l.rstrip().split()
        fields[k] = v
    return fields


def load_matrices(args):
    mats = dd(list)
    for mt in constants.PEAKVIEWER_MATRICES:
        if mt in args:
            mats[mt] = dict(zip(args["%s_names" % mt], [read_mat(fd) for fd in args[mt]]))
    return mats


def render(args, matrices, exons):

    def to_json(value):
        return escape(json.dumps(value))

    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(os.path.join(EXEC_DIR, "../templates/")))
    env.filters.update({'to_json': to_json, 'debug': utils.debug})
    template_file_name = args['view'].replace("-", "_") + "_template.j2"
    sum_template = env.get_template(template_file_name)

    if args['view'] == constants.VIEW_PEAKVIEWER:
        count_pages = 0

        genes = [g.rstrip().split()[0].split(',')[0] for g in args['genes']]

        # Subfolder for summary pages
        summaries_subfolder = "%s/%s" % (args['outdir'], constants.GENES_SUBFOLDER)
        utils.create_if_not_exists(summaries_subfolder)

        while count_pages*constants.MAX_GENES < len(genes):
            prev_page = None
            next_page = None

            subset_genes = genes[count_pages*constants.MAX_GENES: constants.MAX_GENES*(count_pages+1)]
            exons_subset = {k: exons[k] for k in exons.keys() if k in subset_genes}
            name_page = "%03d_%s.html" % (count_pages, args['view'])
            if (count_pages+1)*constants.MAX_GENES < len(genes):
                next_page = "%03d_%s.html" % (count_pages+1, args['view'])
            if not count_pages == 0:
                prev_page = "%03d_%s.html" % (count_pages-1, args['view'])

            full_path = "%s/%s" % (summaries_subfolder, name_page)
            output_file = open(full_path, 'w')
            output_file.write(
                sum_template.render(
                    prevPage=prev_page,
                    nextPage=next_page,
                    namePage=name_page,
                    matrices=matrices,
                    genes=subset_genes,
                    exons=exons_subset

                )
            )
            count_pages += 1

    utils.copyanything(EXEC_DIR+"../static", args['outdir']+"static")
    utils.copyanything(EXEC_DIR+"../GGR-Visual/client/modules", args['outdir']+"static/js/modules/")
    utils.copyanything(EXEC_DIR+"../GGR-Visual/client/stylesheets", args['outdir']+"static/css/")


def main():
    global logger

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
        Time series data peak viewer visualizations
        -------------------------------------------
        '''
    )
    parser.add_argument('-v', action='version', version=constants.VERSION)

    # Base script options
    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument('-o', '--outdir',
                             metavar='outdir',
                             dest='outdir',
                             type=str,
                             required=True,
                             help='Output directory where the files will be placed.')
    base_parser.add_argument('--logger',
                             default=None,
                             help='Path for the logger. Default is output directory')
    base_parser.add_argument('--silent',
                             action='store_true',
                             default=False,
                             help='Silence the logger.')

    # Subparser module to agglutinate all subparsers
    subparsers = parser.add_subparsers(dest='view')
    subparsers.required = True

    # Common options to peak viewers
    common_parser = argparse.ArgumentParser(add_help=False)

    # Peak viewers
    parser_peakviewer = argparse.ArgumentParser(add_help=False)

    parser_peakviewer.add_argument('genes',
                                   type=file,
                                   help='File containing a list of gene IDs to be plotted. ')

    parser_peakviewer.add_argument('exons',
                                   type=file,
                                   help='Column tab-delimited file: <chr>')

    # ChIP-seqs
    parser_peakviewer.add_argument('--tf-matrices', metavar='FACTOR_XX.matrix.tsv [FACTOR_YY.reads.tsv ...]',
                                   type=file,
                                   nargs='*',
                                   help='2 column tab-delimited file: 1st column, gene name; 2nd column reads.')
    parser_peakviewer.add_argument('--tf-matrices-names', metavar='FACTOR_XX [FACTOR_YY ...]',
                                   type=str, nargs='*',
                                   help='ChIP-seq sample names')

    # Chromatin seqs (Histone mods., Chromatin
    parser_peakviewer.add_argument('--histmod-matrices', metavar='HIST_MOD_XX.matrix.tsv [HIST_MOD_YY.matrix.tsv ...]',
                                   type=file,
                                   nargs='*',
                                   help='2 column tab-delimited file: 1st column, gene name; 2nd column reads.')
    parser_peakviewer.add_argument('--histmod-matrices-names', metavar='HIST_MOD_XX [HIST_MOD_YY ...]',
                                   type=str, nargs='*',
                                   help='Histone modification names')

    # # DNaseI-seqs
    # parser_peakviewer.add_argument('--dnase-matrix', metavar='dnaseI.reads.tsv',
    #                                type=file,
    #                                nargs=1,
    #                                help='2 column tab-delimited file: 1st column, gene name; 2nd column reads.')

    subparsers.add_parser(constants.VIEW_PEAKVIEWER,
                          help="Peak viewers for time series data",
                          parents=[base_parser, common_parser, parser_peakviewer])

    # Parse input
    args = parser.parse_args()

    # Create output directory
    # outdir = args.outdir
    if not args.outdir.endswith('/'):
        args.outdir += '/'
    utils.create_if_not_exists(args.outdir)

    # Create logging file
    loggerdir = args.logger
    if loggerdir is None:
        loggerdir = args.outdir
    utils.create_if_not_exists(loggerdir)

    logger = utils.get_logger("%stime-series-viewers.log" % loggerdir, silent=args.silent)
    logger.info("Execution line: %s" % repr(args))

    if args.view == constants.VIEW_PEAKVIEWER:
        exons = read_exons(args.exons)
        matrices = load_matrices(vars(args))
        render(vars(args), matrices, exons)
        
if __name__ == '__main__':
    # Time execution time
    start_time = time.time()
    main()
    end_time = time.time()
    elapsed_str = utils.secs2hms(end_time - start_time)
    logger.info("Execution time: %s" % elapsed_str)
    sys.exit(1)
