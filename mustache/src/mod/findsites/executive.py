import os, sys
from siteFinder import SiteFinder

def action(args, logger):

    genome_path = args['genome']
    bam = args['bam_file']
    peaks_path = args['peak_ranges']

    sitefinder = SiteFinder(genome_path, bam, peaks_path, logger)
    sitefinder.find_sites()