import sys
import shared
from genome import GenomeAlignment
from peakQC import PeakQC

class SiteFinder:

    def __init__(self, genome_path, bam_file, peaks_path, logger=None):

        self.log = logger

        self.bam_file = bam_file
        self.genome_path = genome_path
        self.peaks_path = peaks_path
        self.peaks = shared.get_peaks(peaks_path)

        self.genome = None
        self.peakQC = None

    def find_sites(self):

        self.log_info("Building the genome read depth model...")
        self.genome = GenomeAlignment(self.genome_path, self.bam_file, self.peaks, self.log)
        self.genome.generateModel()
        self.log_debug("Finished building the genome model...")

        self.log_debug("Creating the PeakQC object...")
        self.peakQC = PeakQC(self.genome, self.peaks, self.log)

        for chrom, start_zeroi, end_zeroi in self.peaks:
            peak_str = shared.peak_str(start_zeroi, end_zeroi)

            if self.log: self.log.info("Processing the peak %s..." % peak_str)
            self.peakQC.perform_peak_QC(chrom, start_zeroi, end_zeroi)

        self.peakQC.find_insertion_sites()

        self.peakQC.print_QC_log()


        sys.exit()

            #target_site = self.find_target_site(chrom, start_zeroi, end_zeroi)


    def log_info(self, string):
        if self.log:
            self.log.info(string)


    def log_debug(self, string):
        if self.log:
            self.log.debug(string)