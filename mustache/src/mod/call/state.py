import misc, os
from ..shared.state_superclass import *
from collections import defaultdict

from ..shared.log import Log

class CallState(StateSuperClass):
    def __init__(self, args):
        StateSuperClass.__init__(self)

        self.which = 'call'

        self.paths = CallPaths(fastq_files = args['fastqs'],
                           fastq_to_IS_algnmnts=['fastq1_to_IS.sam', 'fastq2_to_IS.sam'],
                           class_files = args['classifications'],
                           insertion_fasta=args['insertion_fasta'],
                           reference_genomes=args['reference_genomes'],
                           taxon_nodes=args['taxon_nodes'],
                           taxon_names=args['taxon_names'],
                           outdir=args['output_dir'])

        self.settings = CallSettings(args['complete_class_exclusions'],
                                 args['genome_class_exclusions'],
                                 args['insertion_class_exclusions'],
                                 args['threads'],
                                 args['num_reads'])

        self.logger = Log(self.paths.out_dir)

class CallPaths(PathsSuperClass):
    def __init__(self, fastq_files, fastq_to_IS_algnmnts, class_files, insertion_fasta, reference_genomes, taxon_nodes,
                 taxon_names, outdir):
        PathsSuperClass.__init__(self, outdir, taxon_nodes, taxon_names)

        # Directories
        self.sams_dir = self.makedir(os.path.join(outdir, 'sams'))
        self.full_sams_dir = self.makedir(os.path.join(outdir, 'sams', 'full_sams'))
        self.taxon_sorted_sams_dir = self.makedir(os.path.join(outdir, 'sams', 'taxon_sorted_sams_dir'))

        self.peaks_dir = self.makedir(os.path.join(outdir, 'peaks'))
        self.single_peaks_dir = self.makedir(os.path.join(outdir, 'peaks', 'single'))
        self.results_dir = self.makedir(os.path.join(outdir, 'results'))

        # Files
        self.fastq_files = fastq_files
        self.fastq_to_IS_algnmnts = [os.path.join(self.full_sams_dir, os.path.basename(path)) for path in fastq_to_IS_algnmnts]
        self.class_files = class_files
        self.insertion_fasta = insertion_fasta
        self.reference_genomes = reference_genomes
        self.reference_id_to_path = self.create_ref_id_dict(reference_genomes)
        self.fastq_to_genome_algnmnts = None



        self.taxon_sam_paths = {}
        self.taxon_bam_paths = {}
        self.single_peak_paths = {}

        self.indiv_results_out = os.path.join(self.results_dir, 'indiv_results.tsv')
        self.taxonomy_traversal_results_out = os.path.join(self.results_dir, 'taxonomy_traversal_results.tsv')
        self.taxonomy_traversal_subpecies_merged_results_out = os.path.join(self.results_dir, 'taxonomy_traversal_subspecies_merged_results.tsv')

    def makedir(self, path):
        if not os.path.isdir(path):
            os.makedirs(path)
        return os.path.abspath(path)

    def create_ref_id_dict(self, ref_genomes):
        mydict = {}
        for ref in ref_genomes:
            mydict[misc.get_refid_from_path(ref)] = ref

        if len(mydict) != len(ref_genomes):
            raise TypeError("The given reference genomes must have unique base names in their paths.")

        return mydict

class CallSettings(SettingsSuperClass):

    def __init__(self, complete_class_exclusions, genome_class_exclusions, insertion_class_exclusions, threads, num_reads):
        SettingsSuperClass.__init__(self)

        self.complete_class_exclusions = complete_class_exclusions
        self.genome_class_exclusions = genome_class_exclusions
        self.insertion_class_exclusions = insertion_class_exclusions
        self.threads = threads
        self.num_reads = num_reads
        self.taxon_nodes = None
        self.taxon_ranks = None
        self.taxon_names = None
        self.taxon_reads_dict = None
        self.reference_headers_dict = defaultdict(list)
        self.merge_subspecies = True


