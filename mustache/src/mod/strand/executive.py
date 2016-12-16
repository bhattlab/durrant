import os, sys
from strandAnalyzer import StrandAnalyzer

def action(args, logger):

    genome_aln = args['genome_alignment']
    insertion_aln = args['insertion_sequence_alignment']

    genome_rname = args['genome_reference_name']
    range = args['range']
    is_rname = args['insertion_reference_name']

    clip_genome_query_name = args['clip_genome_query_name']


    strand_analyzer = StrandAnalyzer(genome_aln, insertion_aln, genome_rname, range, is_rname, clip_genome_query_name)

    strand_analyzer.run()