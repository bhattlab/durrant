import threading
import time

exitFlag = 0

class filterFlanksThread(threading.Thread):
    def __init__(self, threadID, name, counter):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter
    def run(self):
        print "Starting " + self.name
        print_time(self.name, self.counter, 5)
        print "Exiting " + self.name

def print_time(threadName, delay, counter):
    while counter:
        if exitFlag:
            threadName.exit()
        time.sleep(delay)
        print "%s: %s" % (threadName, time.ctime(time.time()))
        counter -= 1

def filter_flanks_to_fastq(genome_sams, classifs, IS_sams, ref, sorted_dict, state):
    logger = state.logger
    flank_sams_out = {}

    loop_count = 0
    total_read_count = 0
    flanking_reads_count = 0

    for genome_aln1, genome_aln2, IS_aln1, IS_aln2, class1, class2 in izip(genome_sams[0], genome_sams[1],
                                                                           IS_sams[0], IS_sams[1], classifs[0],
                                                                           classifs[1]):
        total_read_count += 1

        loop_count = misc.loop_counter(loop_count, total_read_count, logger)

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
                                             flank_sams_out, sorted_dict, flanking_reads_count)

        elif is_aligned(IS_aln1[0]) and is_aligned(genome_aln2[0]):
            if not passes_genome_class_exclusions(class2, state.settings.genome_class_exclusions):
                continue
            elif not passes_IS_class_exclusions(class1, state.settings.insertion_class_exclusions):
                continue

            flanking_reads_count = write_aln(genome_aln2, IS_aln1, class2, ref, state.paths.taxon_sorted_sams_dir,
                                             flank_sams_out, sorted_dict, flanking_reads_count)


        else:
            continue