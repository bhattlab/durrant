import os


class StateSuperClass:
    def __init__(self):
        pass

class PathsSuperClass:
    def __init__(self, outdir, taxondb):
        self.out_dir = self.makedir(outdir)

        # Files
        self.taxon_nodes = [str(os.path.join(taxondb, 'nodes.dmp')), str(os.path.join(taxondb, 'merged.dmp'))]
        self.taxon_names = str(os.path.join(taxondb, 'names.dmp'))


    def makedir(self, path):
        if not os.path.isdir(path):
            os.makedirs(path)
        return os.path.abspath(path)

class SettingsSuperClass:
    def __init__(self):
        self.path_delim = '-_-mustache-_-'
        self.peak_extension = 1000

        self.nodes_child_parent_overwrite = [('1263037', '47678')]

        self.taxon_names = None
        self.taxon_nodes = None
        self.taxon_ranks = None

    def set_nodes(self, nodes_in):
        for child, parent in self.nodes_child_parent_overwrite:
            nodes_in[child] = parent
        self.taxon_nodes = nodes_in

    def set_names(self, names_in):
        self.taxon_names = names_in

    def set_ranks(self, ranks_in):
        self.taxon_ranks = ranks_in