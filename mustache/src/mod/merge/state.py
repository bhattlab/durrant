import os, json
from os.path import join

from ..shared.log import Log
from ..shared.state_superclass import *

class MergeState(StateSuperClass):
    def __init__(self, args):
        StateSuperClass.__init__(self)

        self.config = dict(json.load(open(args['config'])))
        self.which = 'merge'


        self.paths = MergePaths(mustache_dirs=args['mustache_directories'],
                                output_dir=args['output_dir'],
                                taxondb = self.config['taxondb'])

        self.settings = MergeSettings()

        self.logger = Log(self.paths.out_dir)

class MergePaths(PathsSuperClass):
    def __init__(self, mustache_dirs, output_dir, taxondb):
        PathsSuperClass.__init__(self, output_dir, taxondb)

        # Directories
        self.mustache_dirs = [self.makedir(dir) for dir in mustache_dirs]
        self.merged_bam_dir = self.makedir(join(self.out_dir, 'bams_merged'))
        self.merged_peaks_dir = self.makedir(join(self.out_dir, 'peaks_merged'))
        self.results_dir = self.makedir(join(self.out_dir, 'results'))

        self.bam_info = None
        self.merged_bam_paths = None

        self.merged_peaks_paths = {}

        # Paths
        self.merged_indiv_peaks_path = join(self.results_dir, 'merged_indiv_results.tsv')
        self.merged_taxonomy_traversal = join(self.results_dir, 'merged_taxonomy_traversal_results.tsv')

    def makedir(self, path):
        if not os.path.isdir(path):
            os.makedirs(path)
        return os.path.abspath(path)


class MergeSettings(SettingsSuperClass):
    def __init__(self):
        SettingsSuperClass.__init__(self)