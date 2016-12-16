import multiprocessing
import os
from itertools import izip
from operator import itemgetter
from os.path import basename

import IO
import misc


def sort_flanking_reads(state):

    for ref in state.paths.fastq_to_genome_algnmnts:
        state.logger.info("Filtering IS-flanking reads that mapped to the %s genome..." % basename(ref))
        read_chunks = misc.calc_threads_start_stop(state.settings.num_reads, state.settings.threads)

        process_id = 1
        processes = []
        for start, stop in read_chunks:
            genome_sams = IO.read_genome_alignment(open(state.paths.fastq_to_genome_algnmnts[ref][0]), start), \
                          IO.read_genome_alignment(open(state.paths.fastq_to_genome_algnmnts[ref][1]), start)
            classifs = IO.read_classifications(open(state.paths.class_files[0]), start), \
                       IO.read_classifications(open(state.paths.class_files[1]), start)
            IS_sams = IO.read_insertion_alignments(open(state.paths.fastq_to_IS_algnmnts[0]), start), \
                      IO.read_insertion_alignments(open(state.paths.fastq_to_IS_algnmnts[1]), start)

            filterProcess = filterFlanksProcess(process_id, start, stop, genome_sams,
                                                classifs, IS_sams, ref, state)
            processes.append(filterProcess)
            filterProcess.start()
            process_id += 1
        for proc in processes: proc.join()


class filterFlanksProcess(multiprocessing.Process):

    def __init__(self, processID, start, stop, genome_sams, classifs, IS_sams, ref, state):
        multiprocessing.Process.__init__(self)
        self.processID = processID
        self.startread = start
        self.stopread = stop
        self.genome_sams = genome_sams
        self.classifs = classifs
        self.IS_sams = IS_sams
        self.ref = ref
        self.state = state


    def run(self):
        self.state.logger.debug("Starting filter process number %d..." % self.processID)
        filter_flanks_to_fastq(self.genome_sams, self.classifs, self.IS_sams, self.ref,
                               self.state, self.stopread-self.startread,'PROC%d' % self.processID)



def filter_flanks_to_fastq(genome_sams, classifs, IS_sams, ref, state, maxlines, suffix):
    logger = state.logger
    flank_sams_out = {}

    loop_count = 0
    total_read_count = 0
    flanking_reads_count = 0

    for genome_aln1, genome_aln2, IS_aln1, IS_aln2, class1, class2 in izip(genome_sams[0], genome_sams[1],
                                                                           IS_sams[0], IS_sams[1], classifs[0],
                                                                           classifs[1]):
        if total_read_count == 0:
            state.logger.debug("START READ on %s: %s" % (suffix, class1['HEADER']))

        if total_read_count > maxlines:
            state.logger.debug("STOP READ on %s: %s" % (suffix, class1['HEADER']))
            break
        total_read_count += 1

        loop_count = misc.loop_counter(loop_count, total_read_count, suffix, logger)

        # Check the reads and the classifications align
        check_matching_reads(genome_aln1, genome_aln2, IS_aln1, IS_aln2, class1, class2)
        if not passes_complete_class_exclusions(class1, class2, state.settings.complete_class_exclusions):
            continue

        if is_aligned(IS_aln1[0]) and is_aligned(IS_aln2[0]):
            continue
        elif is_aligned(genome_aln1[0]) and is_aligned(IS_aln2[0]):
            if not passes_genome_class_exclusions(class1, state.settings.genome_class_exclusions):
                continue
            elif not passes_IS_class_exclusions(class2, state.settings.insertion_class_exclusions):
                continue

            flanking_reads_count = write_aln(genome_aln1, IS_aln2, class1, ref, state.paths.taxon_sorted_sams_dir,
                                             flank_sams_out, suffix, flanking_reads_count, state)

        elif is_aligned(IS_aln1[0]) and is_aligned(genome_aln2[0]):
            if not passes_genome_class_exclusions(class2, state.settings.genome_class_exclusions):
                continue
            elif not passes_IS_class_exclusions(class1, state.settings.insertion_class_exclusions):
                continue

            flanking_reads_count = write_aln(genome_aln2, IS_aln1, class2, ref, state.paths.taxon_sorted_sams_dir,
                                             flank_sams_out, suffix, flanking_reads_count, state)

        else:
            continue

    state.logger.info("Total Reads counted on %s: %d" % (suffix, total_read_count))
    state.logger.info("Total Flanking Reads counted on %s: %d" % (suffix, flanking_reads_count))


def write_aln(gen_algnmnts, IS_algnmnts, classif, ref, outdir, sams_out, suffix, flank_count, state):
    refbase = os.path.basename(ref).split('.')[0]
    delim = state.settings.path_delim
    for gen_aln in gen_algnmnts:
        for IS_aln in IS_algnmnts:
            IS = IS_aln['RNAME']
            outfile = os.path.join(outdir, '%s%s%s%s%s%s%s.sam' %
                                   (refbase, delim, classif['TAXID'], delim, IS, delim, suffix))
            if outfile in sams_out:
                outline = [val for key, val in gen_aln.items()]
                outline[0] = outline[0]+':TAXID-%s'%classif['TAXID']
                sams_out[outfile].write("\t".join(outline)+'\n')
                flank_count += 1

            else:
                sams_out[outfile] = open(outfile, 'w')
                outline = [val for key, val in gen_aln.items()]
                outline[0] = outline[0] + ':TAXID-%s' % classif['TAXID']
                sams_out[outfile].write("\t".join(outline) + '\n')
                flank_count += 1
    return flank_count


def sorted_flanks_dict_to_string(sorted_dict, taxon_names=None):
    outlist = []
    for ref in sorted_dict:
        for taxon in sorted_dict[ref]:
            for IS in sorted_dict[ref][taxon]:
                if taxon_names:
                    outlist.append([ref, taxon, taxon_names[taxon], IS, sorted_dict[ref][taxon][IS]])
                else:
                    outlist.append([ref, taxon, IS, sorted_dict[ref][taxon][IS]])
    if taxon_names:
        outlist = sorted(outlist, key=itemgetter(0, 4), reverse=True)
        outlist.insert(0, ['Genome', 'Taxon', 'Name', 'Insertion', 'NFlanks'])
    else:
        outlist = sorted(outlist, key=itemgetter(0, 3), reverse=True)
        outlist.insert(0, ['Genome', 'Taxon', 'Insertion', 'NFlanks'])

    return "\n".join(["\t".join([str(elem) for elem in line]) for line in outlist])

def check_matching_reads(genome_aln1, genome_aln2, IS_aln1, IS_aln2, class1, class2):
    if genome_aln1[0]['QNAME'] == genome_aln2[0]['QNAME'] == IS_aln1[0]['QNAME'] == IS_aln2[0]['QNAME'] == class1['HEADER'] == class2['HEADER']:
        pass
    else:
        raise IndexError("The reads and the classifications need to be in the same order.")

def is_unclassif_read(classif):
    if classif['TAXID'] == '0':
        return True
    return False


def is_aligned(alignment):
    if alignment['RNAME'] == '*':
        return False
    else:
        return True


def passes_complete_class_exclusions(class1, class2, exclusions):
    if exclusions:
        for exclusion in exclusions:
            if class1[exclusion[0]] == exclusion[1] or class2[exclusion[0]] == exclusion[1]:
                return False
    return True

def passes_genome_class_exclusions(classif, exclusions):
    global args
    if exclusions:
        for exclusion in exclusions:
            if classif[exclusion[0]] == exclusion[1]:
                return False
    return True

def passes_IS_class_exclusions(classif, exclusions):
    global args
    if exclusions:
        for exclusion in exclusions:
            if classif[exclusion[0]] == exclusion[1]:
                return False
    return True
