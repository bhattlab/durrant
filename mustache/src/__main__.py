import argparse
import os
from pprint import pformat

from mod.call import executive as call_executive
from mod.merge import executive as merge_executive
from mod.findsites import executive as findsites_executive
from mod.strand import executive as strand_executive


from mod.call import state as call_state
from mod.merge import state as merge_state

from mod.call import argparseTypes as call_argparseTypes
from mod.merge import argparseTypes as merge_argparseTypes
from mod.findsites import argparseTypes as findsites_argparseTypes

from mod.shared.log import SimpleLog

def main(args):

    if args['which'] == 'call':

        callstate = call_state.CallState(args)
        logger = callstate.logger

        logger.info("Here are the arguments as they were given:\n\n%s\n" % pformat(args))

        logger.info("Executing the ISPeaks call protocol...")

        call_executive.action(callstate)

    elif args['which'] == 'merge':
        mergestate = merge_state.MergeState(args)
        logger = mergestate.logger
        logger.info("Here are the arguments as they were given:\n\n%s\n" % pformat(args))

        logger.info("Executing the ISPeaks merge protocol...")

        merge_executive.action(mergestate)

    elif args['which'] == 'findsites':

        logger = SimpleLog()
        logger.info("Here are the arguments as they were given:\n\n%s\n" % pformat(args))
        logger.info("Executing the ISPeaks find-sites protocol...")

        findsites_executive.action(args, logger)

    elif args['which'] == 'strand':

        logger = SimpleLog()
        logger.info("Here are the arguments as they were given:\n\n%s\n" % pformat(args))
        logger.info("Executing the ISPeaks strand protocol...")

        strand_executive.action(args, logger)

    else:
        pass
        #logger.error("Something weird just happened...")



if __name__ == "__main__":

    current_dir = os.path.dirname(__file__)
    data_dir = os.path.join(current_dir, "../data")
    # setup the option parser
    parser = argparse.ArgumentParser(description='')

    subparsers = parser.add_subparsers(help="The argument specifying the type of analysis: single, "
                                            "or merged")

    parser_call = subparsers.add_parser('call', help='Run ISPeaks on a single fastq file to call al lthe peaks')
    parser_call.set_defaults(which="call")

    parser_merge = subparsers.add_parser('merge', help='Run ISPeaks on multiple ISPeaks output folders.')
    parser_merge.set_defaults(which="merge")

    parser_findsite = subparsers.add_parser('find-sites', help='Identify the specific sites of insertion using ISPeaks '
                                                          'insertion-flanking alignments and specific peaks of interest.')
    parser_findsite.set_defaults(which="findsites")

    parser_strand = subparsers.add_parser('strand', help='')
    parser_strand.set_defaults(which="strand")

    # SINGLE arguments
    parser_call.add_argument('-fq', '--fastqs', required=True, nargs = 2, type=call_argparseTypes.fastq_file,
                             help='Two fastq files containing the forward and reverse strands, in that order.')

    parser_call.add_argument('-c', '--classifications', required=True, nargs=2, type=call_argparseTypes.class_file,
                             help='Two tab-separated files where the first column is the read title and the second'
                             'column is the assigned taxon id. The first classification file corresponds to the'
                             'forward reads, and the second file corresponds to the reverse reads, in the same order'
                             'with no reads excluded. MUST HAVE A HEADER.')

    parser_call.add_argument('-cce', '--complete-class-exclusions', required=False, nargs='*',
                             type=call_argparseTypes.complete_class_exclude,
                             help='A filter used to exclude certain reads if either class column matches the value. '
                                    'Input as \"<Column Name1>=<Criterion1> <Column Name2>=<Criterion2>...\"')

    parser_call.add_argument('-gce', '--genome-class-exclusions', required=False, nargs='*',
                             type=call_argparseTypes.genome_class_exclude,
                             default=[call_argparseTypes.genome_class_exclude("PASSEDFILTER=F"),
                                      call_argparseTypes.genome_class_exclude("TAXID=0")],
                             help='A filter used to exclude certain reads if the genome-aligned read has a column '
                                    'that matches the given value. '
                                    'Input as \"<Column Name1>=<Criterion1> <Column Name2>=<Criterion2>...\"')

    parser_call.add_argument('-ice', '--insertion-class-exclusions', required=False, nargs='*',
                             type=call_argparseTypes.IS_class_exclude,
                             help='A filter used to exclude certain reads if the IS-aligned read has a column '
                                    'that matches the given value. '
                                    'Input as \"<Column Name1>=<Criterion1> <Column Name2>=<Criterion2>...\"')

    parser_call.add_argument('-o', '--output-dir', required=True, type=call_argparseTypes.output_folder,
                             help='Specify the output folder to create. If already created, it must be empty.')

    parser_call.add_argument('-is', '--insertion-fasta', required=True, type=call_argparseTypes.insertion_fasta,
                             default=os.path.join(data_dir, "IS_fastas/Bacteroides_all.fasta"),
                             help='A fasta file containing the insertion sequences of interest,'
                             ' concatenated sequentially in any order.')

    parser_call.add_argument('-g', '--reference-genomes', nargs='+',
                             type=call_argparseTypes.genome, required = True,
                             help='Reference genomes to analyze in fasta format. Input must be '
                                'in the format \"<reference-fasta1> <reference-fasta2>\"')

    parser_call.add_argument('-nodes', '--taxon-nodes', required=False, type=call_argparseTypes.taxon_nodes,
                             default=[call_argparseTypes.taxon_nodes(os.path.abspath(os.path.join(data_dir, "TaxonomyDatabase/nodes.dmp"))),
                                      call_argparseTypes.taxon_nodes(os.path.abspath(os.path.join(data_dir, "TaxonomyDatabase/merged.dmp")))],
                             help='Location of the NCBI Taxonomy Database nodes.dmp and/or merged.dmp files',
                             nargs='+')

    parser_call.add_argument('-names', '--taxon-names', required=False, type=call_argparseTypes.taxon_names,
                             default=call_argparseTypes.taxon_nodes(
                                   os.path.abspath(os.path.join(data_dir, "TaxonomyDatabase/names.dmp"))),
                             help='Location of the NCBI Taxonomy Database names.dmp')

    parser_call.add_argument('-n', '--num_reads', required=False, type=int,
                               help='The total number of paired-end reads. If not given, it will be calculated by counting'
                                    'the number of lines in one of the given classification files')

    parser_call.add_argument('-p', '--threads', required=False,
                        default=1, type = int,
                        help='The number of threads to run with bowtie2 alignments. This is typically (number'
                             ' of available processors) / 2.')

    # MERGED arguments
    parser_merge.add_argument('-d', '--ispeaks_directories', required=True,
                              nargs='+', action=merge_argparseTypes.minimum_length(2),
                              type=merge_argparseTypes.ispeaks_directory,
                              help='Include 2 or more ISPeaks output directories to merge')

    parser_merge.add_argument('-o', '--output-dir', required=True, type=merge_argparseTypes.output_folder,
                              help='Specify the output folder to create. If already created, it must be empty.')

    parser_merge.add_argument('-nodes', '--taxon-nodes', required=False, type=call_argparseTypes.taxon_nodes,
                             default=[call_argparseTypes.taxon_nodes(
                                 os.path.abspath(os.path.join(data_dir, "TaxonomyDatabase/nodes.dmp"))),
                                      call_argparseTypes.taxon_nodes(
                                          os.path.abspath(os.path.join(data_dir, "TaxonomyDatabase/merged.dmp")))],
                             help='Location of the NCBI Taxonomy Database nodes.dmp and/or merged.dmp files',
                             nargs='+')

    parser_merge.add_argument('-names', '--taxon-names', required=False, type=call_argparseTypes.taxon_names,
                             default=call_argparseTypes.taxon_nodes(
                                 os.path.abspath(os.path.join(data_dir, "TaxonomyDatabase/names.dmp"))),
                             help='Location of the NCBI Taxonomy Database names.dmp')

    # Find-Site arguments
    parser_findsite.add_argument('-b', '--bam-file', required=True,
                              type=findsites_argparseTypes.bam_file,
                              help='The IS-flanking alignment file in bam format.')
    parser_findsite.add_argument('-p', '--peak-ranges', required=True,
                              type=findsites_argparseTypes.peak_ranges,
                              help='A tab-delimited text file containing three columns of genomic positions corresponding'
                                   'to peak ranges of interest.'
                                   'The left column includes the specific chromosome reference, the middle position must precede the right column positions.')
    parser_findsite.add_argument('-g', '--genome', required=True,
                                 type=findsites_argparseTypes.genome)


    ## Strand-Orient
    parser_strand.add_argument('-ga', '--genome-alignment', required=True)
    parser_strand.add_argument('-isa', '--insertion-sequence-alignment', required=True)
    parser_strand.add_argument('-gsr', '--genome-reference-name', required=True, type=str)
    parser_strand.add_argument('-r', '--range', required=True, nargs=2, type=int)
    parser_strand.add_argument('-isr', '--insertion-reference-name', required=True, type=str)
    parser_strand.add_argument('--clip-genome-query-name', required=False, default=True, type = bool)






    args = parser.parse_args()
    args = vars(args)

    main(args)