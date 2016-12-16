class StateSuperClass:
    def __init__(self):
        pass

class PathsSuperClass:
    def __init__(self, outdir, taxon_nodes, taxon_names):
        self.out_dir = self.makedir(outdir)

        # Files
        self.taxon_nodes = taxon_nodes
        self.taxon_names = taxon_names

class SettingsSuperClass:
    def __init__(self):
        self.path_delim = '-_-ispeaks-_-'
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