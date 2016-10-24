import os
import sys
import argparse
import time
import re
import h5py
import simplejson
import utils
import constants
from collections import defaultdict as dd
from jinja2 import Environment, FileSystemLoader, escape
import json

EXEC_DIR = utils.module_path() + "/"


def read_exons(exonsf):
    res = {}
    for l in exonsf:
        if l.startswith('#'): continue
        chrom, strand, gene, exons = l.rstrip().split('\t')
        gene_name = None
        if '(' in gene:  # TODO: temporary solution for current files -- REMOVE
            gene_name = gene.split('(')[0]
            gene = re.match(".*\((.*)\).*", gene).groups()[0]
        res[gene] = {'chrom': chrom, 'strand': strand, 'exons': exons, 'gene_name': gene_name}
    return res


def read_mat(mat):
    return h5py.File(mat.name, 'r')


def load_matrices(args):
    """Load data from hdf5 files"""
    mats = dd(list)
    res = []
    for mt in constants.PEAKVIEWER_MATRICES:
        if args[mt]:
            factor_name = args["%s_names" % mt]
            if not factor_name:
                factor_name = [os.path.basename(factor_path.name).split('.')[0] for factor_path in args[mt]]
            mats[mt] = dict(zip(factor_name, [read_mat(fd) for fd in args[mt]]))
            res_factor = [list(mm[1]['resolutions']) for mm in mats[mt].items()]
            if res and set(res) != set(res_factor):
                logger.error("")
            res.extend(res_factor)
    return mats, res[0]


def render(args, matrices, exons, resolution, timepoints):

    def to_json(value):
        return escape(json.dumps(value))

    def nparray_to_json(value):
        return escape(simplejson.dumps(value.tolist()))

    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(os.path.join(EXEC_DIR, "../templates/")))
    env.filters.update({'to_json': to_json, 'nparray_to_json': nparray_to_json, 'debug': utils.debug})
    template_file_name = args['view'].replace("-", "_") + "_template.j2"
    sum_template = env.get_template(template_file_name)

    if args['view'] == constants.VIEW_PEAKVIEWER:
        count_pages = 0

        genes = exons.keys()

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
                    exons=exons_subset,
                    resolutions=resolution,
                    timepoints=timepoints
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

    # parser_peakviewer.add_argument('genes',
    #                                type=file,
    #                                help='File containing a list of gene IDs to be plotted. ')

    parser_peakviewer.add_argument('exons',
                                   type=file,
                                   help='Column tab-delimited file: <chrom> <strand> <gene> <exons>')

    # TODO: extract the timepoints from HDF5 files (first I have to add them.., this is it for now)
    parser_peakviewer.add_argument('--timepoints',
                                   type=float, nargs='+', default=[0],
                                   help='Time points in hours (0 --> 0hs, 0.5 --> 30min, etc.')

    # ChIP-seqs
    parser_peakviewer.add_argument('--tfs', metavar='FACTOR_XX.h5 [FACTOR_YY.h5...]',
                                   type=file,
                                   nargs='*',
                                   help='2 column tab-delimited file: 1st column, gene name; 2nd column reads.')

    # Chromatin seqs (Histone mods., Chromatin
    parser_peakviewer.add_argument('--histmods', metavar='HIST_MOD_XX.h5 [HIST_MOD_YY.h5 ...]',
                                   type=file,
                                   nargs='*',
                                   help='2 column tab-delimited file: 1st column, gene name; 2nd column reads.')

    # # DNaseI-seqs
    parser_peakviewer.add_argument('--dnase-matrix', metavar='dnaseI.reads.tsv',
                                   type=file,
                                   nargs=1,
                                   help='2 column tab-delimited file: 1st column, gene name; 2nd column reads.')

    subparsers.add_parser(constants.VIEW_PEAKVIEWER,
                          help="Peak viewers for time series data",
                          parents=[base_parser, common_parser, parser_peakviewer])

    # Parse input
    args = parser.parse_args()

    if not args.tfs and not args.histmods:
        print "ERROR :: No matrix input specified"
        sys.exit(1)

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
        matrices = None
        try:
            matrices, res = load_matrices(vars(args))
            render(vars(args), matrices, exons, res, args.timepoints)
        finally:
            if matrices:
                [mm.close() for m in matrices.values() for mm in m.values()]


if __name__ == '__main__':
    # Time execution time
    start_time = time.time()
    main()
    end_time = time.time()
    elapsed_str = utils.secs2hms(end_time - start_time)
    logger.info("Execution time: %s" % elapsed_str)
    sys.exit(1)
