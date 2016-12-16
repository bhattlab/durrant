import pysam
import numpy
import sys
import operator
import shared
from collections import defaultdict
from Bio import SeqIO

class GenomeAlignment:

    def __init__(self, genome_path, bam_path, peaks, logger=None):

        self.log = logger

        self.bam_path = bam_path
        self.peaks = peaks
        self.genome_seq = self.get_genome(genome_path)

        self.chrom_lengths = self.get_chrom_lengths(bam_path)

        self.total_depth = {}
        self.valid_depth = {}

        self.forward_depth = {}
        self.reverse_depth = {}
        self.evenly_mixed_depth = {}

        self.forward_clipped_depth = {}
        self.reverse_clipped_depth = {}


    def get_chrom_lengths(self, bam_path):

        chrom_lengths = {}
        for SQ in pysam.AlignmentFile(bam_path, 'rb').header['SQ']:
            chrom_lengths[SQ['SN']] = SQ['LN']
        if len(chrom_lengths) == 0:
            raise TypeError("The header file does not contain the appropriate chromosome lengths.")

        return chrom_lengths


    def generateModel(self):
        self.log_info("Building the chromosome position arrays in preparation for read processing...")
        self.build_chromosome_position_lists()

        self.log_info("Processing the reads in the given alignment file...")
        self.process_reads()
        self.log_info("%d bases contain at least 1 read..." % self.get_total_aligned_bp(self.total_depth))
        self.log_info("%d bases contain at least 1 valid read..." % self.get_total_aligned_bp(self.valid_depth))
        self.log_info("%d bases contain at least 1 forward-oriented read..." % self.get_total_aligned_bp(self.forward_depth))
        self.log_info("%d bases contain at least 1 reverse-oriented read..." % self.get_total_aligned_bp(self.reverse_depth))
        self.log_info("%d bases contain relatively even proportions of forward- and reverse-oriented reads" %
                      self.get_total_aligned_bp(self.evenly_mixed_depth))
        self.log_info("%d bases contain at least 1 forward-oriented soft clipping..." % self.get_total_aligned_bp(self.forward_clipped_depth))
        self.log_info("%d bases contain at least 1 reverse-oriented soft_clipping..." % self.get_total_aligned_bp(self.reverse_clipped_depth))



    def build_chromosome_position_lists(self):

        for chrom in self.chrom_lengths:
            self.total_depth[chrom] = [0]*self.chrom_lengths[chrom]
            self.valid_depth[chrom] = [0] * self.chrom_lengths[chrom]

            self.forward_depth[chrom] = [0]*self.chrom_lengths[chrom]
            self.reverse_depth[chrom] = [0] * self.chrom_lengths[chrom]
            self.evenly_mixed_depth[chrom] = [0] * self.chrom_lengths[chrom]

            self.forward_clipped_depth[chrom] = [0] * self.chrom_lengths[chrom]
            self.reverse_clipped_depth[chrom] = [0] * self.chrom_lengths[chrom]


    def process_reads(self):

        for read in pysam.AlignmentFile(self.bam_path, 'rb'):
            true_start = read.reference_start-read.query_alignment_start
            true_end = true_start + read.query_length
            chrom = read.reference_name

            for i in range(true_start, true_end):
                try:

                    self.total_depth[chrom][i] += 1
                    if not read.is_reverse:
                        if i >= read.reference_start and i < read.reference_end:
                            self.forward_depth[chrom][i] += 1
                            self.valid_depth[chrom][i] += 1

                        else:
                            self.forward_clipped_depth[chrom][i] += 1

                    else:
                        if i >= read.reference_start and i < read.reference_end:
                            self.reverse_depth[chrom][i] += 1
                            self.valid_depth[chrom][i] += 1
                        else:
                            self.reverse_clipped_depth[chrom][i] += 1


                except IndexError:
                    continue

        ### Now focus on peaks to count stacked sites
        for chrom, start, end in self.peaks:
            for pos in range(start,end+1):
                try:
                    ratio = float(self.forward_depth[chrom][pos]) / float(self.reverse_depth[chrom][pos])
                except ZeroDivisionError:
                    continue
                if 0.1 <= ratio <= 1.9:
                    both_count = self.forward_depth[chrom][pos] + self.reverse_depth[chrom][pos]
                    self.evenly_mixed_depth[chrom][pos] = both_count


    def print_read_counts(self):
        print '\t'.join(["CHROM", "POS", "TOTAL_DEPTH","VALID_DEPTH",
                         "FORWARD_DEPTH", "REVERSE_DEPTH", "EVENLY_MIXED_DEPTH",
                         "FORWARD_CLIPPED_DEPTH", "REVERSE_CLIPPED_DEPTH"])

        for chrom in self.chrom_lengths:
            for pos in range(self.chrom_lengths[chrom]):
                print '\t'.join([chrom, str(pos+1),
                                 str(self.total_depth[chrom][pos]),
                                 str(self.valid_depth[chrom][pos]),
                                 str(self.forward_depth[chrom][pos]),
                                 str(self.reverse_depth[chrom][pos]),
                                 str(self.evenly_mixed_depth[chrom][pos]),
                                 str(self.forward_clipped_depth[chrom][pos]),
                                 str(self.reverse_clipped_depth[chrom][pos])])

    ## Simple Calculations
    def get_total_aligned_bp(self, depth_dict):
        total_aligned_bp = 0
        for chrom in depth_dict:
            chrom_length = len(depth_dict[chrom])
            aligned_bp = chrom_length - depth_dict[chrom].count(0)
            total_aligned_bp += aligned_bp
        return total_aligned_bp

    def get_genome(self, genome_path):
        genome_seq = {}
        for record in SeqIO.parse(genome_path, 'fasta'):
            genome_seq[record.id] = str(record.seq)
        return genome_seq

    ## Logging Wrappers
    def log_info(self, string):
        if self.log:
            self.log.info(string)

    def log_debug(self, string):
        if self.log:
            self.log.debug(string)









