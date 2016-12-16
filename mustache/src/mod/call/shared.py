import sys

def get_taxon_nodes_ranks(nodes_locations, logger=None):
    assert type(nodes_locations) is list, "The nodes location must be a list of file locations."

    taxon_nodes_dict = {}
    taxon_ranks_dict = {}
    for location in nodes_locations:
        with open(location) as nodes_in:
            for line in nodes_in:
                line = line.strip().split("|")
                id = line[0].strip()
                rank = line[2].strip()
                parent_id = line[1].strip()
                taxon_nodes_dict[id] = parent_id
                taxon_ranks_dict[id] = rank
    return taxon_nodes_dict, taxon_ranks_dict

def get_taxon_names(names_location, logger=None):
    assert type(names_location) is str, "The nodes location must be a string to the file location."

    taxon_names_dict = {}
    with open(names_location) as names_in:
        for line in names_in:
            line = line.strip().split("|")
            name_type = line[3].strip()
            taxon_id = line[0].strip()
            name = line[1].strip()
            if name_type == 'scientific name':
                taxon_names_dict[taxon_id] = name

    return taxon_names_dict

def get_taxon_hierarchy_list(taxon_id, taxon_nodes_dict):
    hierarchy = [taxon_id]

    while taxon_id != '1' and taxon_id != '0':
        taxon_id = taxon_nodes_dict[taxon_id]
        hierarchy.append(taxon_id)

    return hierarchy

def get_taxon_hierarchy_set(taxon_id, taxon_nodes_dict):
    hierarchy = set([taxon_id])

    while taxon_id != '1' and taxon_id != '0':
        taxon_id = taxon_nodes_dict[taxon_id]
        hierarchy.add(taxon_id)

    return hierarchy

def is_parent_child(start_taxon1, hierarchy1, start_taxon2, hierarchy2):
    if (start_taxon1 not in hierarchy2) and (start_taxon2 in hierarchy1):
        return True
    if (start_taxon2 not in hierarchy1) and (start_taxon1 in hierarchy2):
        return True
    return False

def which_parent_child(start_taxon1, hierarchy1, start_taxon2, hierarchy2):
    if (start_taxon1 not in hierarchy2) and (start_taxon2 in hierarchy1):
        return 0
    if (start_taxon2 not in hierarchy1) and (start_taxon1 in hierarchy2):
        return 1
    return False


