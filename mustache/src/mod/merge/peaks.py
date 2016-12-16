import os
import sys
from collections import defaultdict
from glob import glob
from operator import itemgetter
from os.path import basename
from pprint import pprint

import pysam

import IO
import macs2
import misc


def call_all_peaks(version, bamdict, peak_paths, outdir, state):
    for key in bamdict:
        state.logger.debug("CALLING PEAKS FOR %s" % key)
        bam = bamdict[key]
        handle = pysam.AlignmentFile(bam, "rb")
        genome_len = handle.header['SQ'][0]['LN']
        read_count = sum([1 for read in handle])
        state.logger.debug("READ COUNT %s" % read_count)
        if read_count < 10:
            state.logger.info("Too few reads to perform peak calling in %s" %
                              "-".join(basename(bam).split(state.settings.path_delim)))
            continue
        state.logger.debug("PASSED READ THRESHOLD")
        peak_path = macs2.call_peaks(version, bam, genome_len, outdir)
        state.logger.debug("CALLED THE PEAKS")
        if check_peak_path(peak_path):
            peak_paths[key] = peak_path
        else:
            delete_peak_files(peak_path)

def check_peak_path(path):
    if len(open(path).read()) > 0:
        return True
    else:
        return False

def delete_peak_files(path):
    prefix = ".".join(path.split('.')[:-1])
    for peakfile in glob(prefix+'*'):
        os.remove(peakfile)


def process_peaks_indiv(merged_peak_paths, orig_sam_info, outfile, state):
    final_results = []

    # Getting the peaks for the starting nodes

    for key in merged_peak_paths:
        read_count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        genome, taxon, IS = key.split(state.settings.path_delim)
        path = merged_peak_paths[key]

        peaks = get_peaks(path)

        merged_peaks = {}
        for chrom in peaks:
            merged_peaks[chrom] = merge_peaks_single(peaks[chrom], state.settings.peak_extension)

        for chrom in merged_peaks:
            #print chrom
            for peak in merged_peaks[chrom]:
                #print '\t',peak
                for sample in orig_sam_info[key]:
                    #print '\t','\t',sample
                    for read in pysam.AlignmentFile(orig_sam_info[key][sample], 'rb'):
                        if is_within_peak(read.pos, peak[0], peak[1], state.settings.peak_extension):
                            read_count_dict[sample][chrom]['-'.join([str(i) for i in peak])] += 1

        for sample in read_count_dict:
            for chrom in read_count_dict[sample]:
                for peak in read_count_dict[sample][chrom]:
                    final_results.append([sample, genome, taxon, state.settings.taxon_names[taxon], IS, chrom,
                                          peak.split('-')[0], peak.split('-')[1], str(read_count_dict[sample][chrom][peak])])

    final_results.insert(0, ['Sample', 'Genome', 'TaxonID', 'TaxonName', 'Insertion', 'Chrom', 'PeakStart', 'PeakEnd',
                             'FlankingReadCount'])
    open(outfile, 'w').write('\n'.join(['\t'.join(line) for line in final_results]))


def process_peaks_taxonomy_traversal(merged_peak_paths, orig_bam_info, outfile, state):
    final_results = []

    # Getting the peaks for the starting nodes
    for key1 in merged_peak_paths:
        genome1, taxon1, IS1 = key1.split(state.settings.path_delim)
        taxonomy = misc.get_taxon_hierarchy_set(taxon1, state.settings.taxon_nodes)

        parent_peaks_keys = get_parent_peaks_keys(genome1, taxon1, IS1, taxonomy, merged_peak_paths, state)

        if len(parent_peaks_keys) > 0:
            merged_peaks = merge_peaks_taxonomy(genome1, taxon1, IS1, parent_peaks_keys, merged_peak_paths, state)

            samples = get_samples(orig_bam_info)
            peak_dict = create_peak_dict(merged_peaks, samples)

            # Now assign the reads in the sam files to specified read bins...
            peak_dict = assign_reads_to_peaks(orig_bam_info[key1], merged_peaks, peak_dict, state.settings.peak_extension)

            for key in parent_peaks_keys:
                peak_dict = assign_reads_to_peaks(orig_bam_info[key], merged_peaks, peak_dict,
                                                  state.settings.peak_extension)

            # Write the peak counts to the final results
            parent_taxa = [key.split(state.settings.path_delim)[1] for key in parent_peaks_keys]
            out_taxa = [taxon1] + parent_taxa
            for sample in samples:
                for chrom in peak_dict[sample]:
                    for peak in peak_dict[sample][chrom]:
                        final_results.append([sample,
                                              genome1,
                                              taxon1,
                                              state.settings.taxon_names[taxon1].replace(' ','_'),
                                              IS1, chrom, peak.split('-')[0], peak.split('-')[1],
                                              str(peak_dict[sample][chrom][peak])])

        # SINGLE NODE - NO TRAVERSAL NEEDED - WAS COMPLETED ALREADY IN
        #  process_peaks_indiv()
        else:
            continue

    final_results.insert(0, ['Sample', 'Genome', 'StartTaxonIDs', 'StartTaxonName', 'Insertion', 'Chrom', 'PeakStart', 'PeakEnd', 'FlankingReadCount'])
    open(outfile, 'w').write('\n'.join(['\t'.join(line) for line in final_results]))


def assign_reads_to_peaks(bam_info, merged_peaks, peak_dict, extension=0):

    for sample in bam_info:
        for read in pysam.AlignmentFile(bam_info[sample], 'rb'):

            for chrom in merged_peaks:

                for peak in merged_peaks[chrom]:

                    if is_within_peak(read.pos, int(peak[0]), int(peak[1]), extension):
                        peak_dict[sample][chrom]['-'.join([str(i) for i in peak])] += 1

    return peak_dict

def create_peak_dict(peaks_in, samples):
    peak_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for sample in samples:
        for chrom in peaks_in:
            for peak in peaks_in[chrom]:
                peak_dict[sample][chrom]["-".join([str(i) for i in peak])] = 0

    return peak_dict

def get_samples(bam_info):
    samples = set()
    for key in bam_info:
        samples = samples | set(bam_info[key].keys())
    return samples


def get_parent_peaks_keys(genome1, taxon1, IS1, taxonomy, peak_paths, state):
    parent_peaks = []
    for key2 in peak_paths:
        genome2, taxon2, IS2 = key2.split(state.settings.path_delim)
        if genome1 == genome2 and IS1 == IS2 and taxon2 in taxonomy:
            parent_peaks.append(key2)
    return parent_peaks


def merge_peaks_taxonomy(genome1, taxon1, IS1, parent_peaks_keys, peak_paths, state):
    child_peaks_path = peak_paths[state.settings.path_delim.join([genome1, taxon1, IS1])]
    child_peaks = get_peaks(child_peaks_path)
    #print "CHILD PEAKS:", genome1, taxon1, IS1
    #pprint(child_peaks)

    merged_peaks = child_peaks
    for key in parent_peaks_keys:
        parent_genome, parent_taxon, parent_IS = key.split(state.settings.path_delim)
        parent_peaks_path = peak_paths[state.settings.path_delim.join([parent_genome, parent_taxon, parent_IS])]
        parent_peaks = get_peaks(parent_peaks_path)
        #print "PARENT PEAKS:", parent_genome, parent_taxon, parent_IS
        #pprint(parent_peaks)

        tmp_peaks = merged_peaks
        for chrom in merged_peaks:
            try:
                peaks1, peaks2 = merged_peaks[chrom], parent_peaks[chrom]
                tmp_peaks[chrom] = merge_peaks_double_left(peaks1, peaks2, state.settings.peak_extension)
            except KeyError:
                continue
        merged_peaks = tmp_peaks

    #print "MERGED PEAKS:"
    #pprint(merged_peaks)
    return merged_peaks

def get_peaks(peaks_path):
    peaks = defaultdict(list)
    with open(peaks_path) as infile:
        for line in infile:
            line = line.strip().split()
            genome_frag, start, stop = line[0], int(line[1]), int(line[2])
            peaks[genome_frag].append([start, stop])
    return dict(peaks)


def merge_peaks_single(peaks, extension=0):
    if len(peaks) == 0 or len(peaks) == 1: return peaks
    peaks = peak_gen(peaks)


    peaks_merged = []
    peak = peaks.next()
    peak = [peak[0] - extension, peak[1] + extension]

    peak_next = [-sys.maxint, -sys.maxint]
    while True:
        try:
            peak_next = peaks.next()
            peak_next = [peak_next[0] - extension, peak_next[1] + extension]
            #print peak_next
            overlap = False
            while is_overlap(peak, peak_next):
                overlap = True
                try:
                    peak = [min([peak_next[0], peak[0]]), max([peak_next[1], peak[1]])]
                    peak_next = peaks.next()
                    peak_next = [peak_next[0] - extension, peak_next[1] + extension]

                except StopIteration:
                    break
            peaks_merged.append([peak[0]+extension, peak[1]-extension])
            peak = peak_next


        except StopIteration:
            break

    if not is_overlap(peaks_merged[-1], [peak_next[0]+extension, peak_next[1]-extension]):
        peaks_merged.append([peak_next[0] + extension, peak_next[1] - extension])

    return peaks_merged

def merge_peaks_double_left(peaks1, peaks2, extension=0):
    peaks_double = sorted(peaks1 + peaks2, key=itemgetter(0))
    peaks_merged = merge_peaks_single(peaks_double, extension)

    peaks_checked = []
    for peak1 in peaks1:
        for peak2 in peaks_merged:
            if is_overlap(peak1, peak2, extension):
                peaks_checked.append(peak2)

    final_peaks = merge_peaks_single(peaks_checked, extension)

    return final_peaks

def peak_gen(peaks):
    for peak in peaks:
        yield peak

def is_overlap(peak1, peak2, extension=0):
    if peak1[0]-extension <= peak2[1]+extension and peak2[0]-extension <= peak1[1]+extension:
        return True
    else:
        return False


def is_within_peak(pos, peak_start, peak_end, extension):
    dist_to_peak = sys.maxint

    if pos > peak_start - extension and pos < peak_end + extension:
        return True

    return False
