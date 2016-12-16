import os, sys
import time
from os.path import basename
from collections import defaultdict
from glob import glob
import samtools, peaks
import shared

def action(state):

    logger = state.logger
    logger.info("Loading taxonomy nodes and ranks from specified file...")
    nodes, ranks = shared.get_taxon_nodes_ranks(state.paths.taxon_nodes)
    state.settings.set_nodes(nodes)
    state.settings.set_ranks(ranks)

    logger.info("Loading taxonomy names from specified file...")
    state.settings.set_names(shared.get_taxon_names(state.paths.taxon_names))

    logger.info("Collecting all bam file information...")
    collect_bam_info(state)

    logger.info("Merging all of the concordant bam files between samples...")
    samtools.merge_bam_files(state)

    logger.info("Beginning peak calling...")
    logger.info("Sorting and indexing the bam files...")
    samtools.process_all_taxon_sams(state)

    logger.info("Calling insertion peaks for each individual merged bam file...")
    peaks.call_all_peaks('2.1.1.20160309', state.paths.merged_bam_paths, state.paths.merged_peaks_paths,
                         state.paths.merged_peaks_dir, state)

    #logger.info("Processing the merged peaks at an individual taxon level...")
    #peaks.process_peaks_indiv(state.paths.merged_peaks_paths, state.paths.bam_info, state.paths.merged_indiv_peaks_path,
    #                          state)
    #logger.info("Saved the results to %s" % basename(state.paths.merged_indiv_peaks_path))


    logger.info("Processing the merged peaks by taxonomy traversal binning...")
    peaks.process_peaks_taxonomy_traversal(state.paths.merged_peaks_paths, state.paths.bam_info,
                                           state.paths.merged_taxonomy_traversal, state)
    logger.info("Saved the results to %s" % basename(state.paths.merged_taxonomy_traversal))

    logger.info("Analysis Complete :)")

def collect_bam_info(state):

    bam_info = defaultdict(lambda: defaultdict(str))
    for dir in state.paths.ispeaks_dirs:
        sam_dir = os.path.join(dir, 'sams', 'taxon_sorted_sams_dir')
        if not os.path.isdir(sam_dir):
            raise TypeError("The directory architecture of the ISPeaks folder %s is unexpected." % basename(dir))

        for sam in glob(sam_dir+'/*.bam'):
            key = '.'.join(basename(sam).split('.')[:-1])
            bam_info[key][dir] = sam

    state.paths.bam_info = bam_info

