import os
from collections import defaultdict
from glob import glob
from os.path import basename

from Bio import SeqIO

import IO


def count_lines(file_in, skip_header=False):
    count = 0
    for line in file_in: count += 1

    if skip_header: return count-1
    else: return count

def loop_counter(loop_count, total_read_count, suffix, logger, loop_count_out=1000000):
    loop_count += 1
    if loop_count == loop_count_out:
        logger.info("Total Reads Processed on %s: %d" % (suffix, total_read_count))
        loop_count = 0
    return loop_count

def calc_threads_start_stop(num_reads, threads):
    read_pos = 1
    reads_per_thread = num_reads / threads
    remainder = num_reads-(threads*reads_per_thread)

    start_stop = []
    for i in range(threads):
        if i != threads - 1:
            start_stop.append((read_pos, read_pos+reads_per_thread-1))
        else:
            start_stop.append((read_pos, read_pos + reads_per_thread - 1 + remainder))
        read_pos += reads_per_thread

    return start_stop

def consolidate_taxon_sams(state):
    samsdir = state.paths.taxon_sorted_sams_dir
    delim = state.settings.path_delim
    sams = glob(samsdir+'/*')
    sams_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for sam in sams:
        genome, taxon, insert, proc = basename(sam).split(state.settings.path_delim)
        readcount = count_lines(open(sam), skip_header = False)
        sams_dict[genome][taxon][insert] += readcount

    state.settings.taxon_reads_dict = sams_dict
    for genome in sams_dict:

        header = state.settings.reference_headers_dict[genome][0]

        for taxon in sams_dict[genome]:

            for IS in sams_dict[genome][taxon]:
                match = os.path.join(samsdir, delim.join([genome, taxon, IS]) + delim)
                group = glob(match+'*')
                outfile = match.strip(delim)+'.sam'
                cat_delete(outfile, group, header)
                state.paths.taxon_sam_paths[state.settings.path_delim.join([genome, taxon, IS])] = outfile

def cat_delete(outpath, group, header=None):

    with open(outpath,'w') as out:
        if header:
            out.write('\n'.join(header)+'\n')
        for path in group:
            with open(path) as infile:
                for line in infile:
                    out.write(line)
            os.remove(path)

def get_genome_length(genome_path):
    length = 0
    with open(genome_path, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            length += len(record.seq)
    return length

def get_refid_from_path(genome_path):
    return basename(genome_path).split('.')[0]

def get_sam_header(path):
    with open(path) as infile:
        header = []
        line = infile.readline().strip()
        while line[0] == '@':
            header.append(line)
            line = infile.readline().strip()
    return header

def get_taxon_hierarchy_set(taxon_id, taxon_nodes_dict, exclude_orig=True):
    if exclude_orig == True:
        hierarchy = set()
    else:
        hierarchy = set([taxon_id])

    while taxon_id != '1' and taxon_id != '0':
        taxon_id = taxon_nodes_dict[taxon_id]
        hierarchy.add(taxon_id)

    return hierarchy


def merge_subspecies_sams(state):
    merge_dict = {}
    for key in state.paths.taxon_sam_paths:
        genome, taxon, IS = key.split(state.settings.path_delim)
        if state.settings.taxon_ranks[taxon] == 'species':
            merge_dict.update({key:[]})

    for key1 in merge_dict:
        genome1, taxon1, IS1 = key1.split(state.settings.path_delim)
        for key2 in state.paths.taxon_sam_paths:
            genome2, taxon2, IS2 = key2.split(state.settings.path_delim)
            if genome1 == genome2 and IS1 == IS2 and is_subspecies(taxon1, taxon2, state.settings.taxon_nodes):
                merge_dict[key1].append(key2)

    merged_sam_paths = set()
    remove_subspecies_keys = set()
    for key in merge_dict:
        genome, taxon, IS = key.split(state.settings.path_delim)
        if len(merge_dict[key]) > 0:

            merged_sam_path = state.paths.taxon_sam_paths[key]+'.tmp'
            merged_sam_paths.add(merged_sam_path)

            with open(merged_sam_path, 'w') as out:
                out.write('\n'.join(state.settings.reference_headers_dict[genome][0])+'\n')
                for read in IO.read_genome_alignment(open(state.paths.taxon_sam_paths[key]), 1):
                    for aln in read:
                        out.write('\t'.join(aln.values())+'\n')
                for subspecies in merge_dict[key]:
                    remove_subspecies_keys.add(subspecies)
                    for read in IO.read_genome_alignment(open(state.paths.taxon_sam_paths[subspecies]), 1):
                        for aln in read:
                            out.write('\t'.join(aln.values()) + '\n')

    for key in remove_subspecies_keys:
        os.remove(state.paths.taxon_sam_paths[key])
        state.paths.taxon_sam_paths.pop(key, None)

    for path in merged_sam_paths:
        old = '.'.join(path.split('.')[:-1])
        os.rename(path, old)

def is_subspecies(species, query, taxon_nodes):
    taxonomy = get_taxon_hierarchy_set(query, taxon_nodes)
    if species in taxonomy:
        return True
    else:
        return False