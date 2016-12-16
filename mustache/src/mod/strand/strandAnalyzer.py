import pysam

class StrandAnalyzer:
    def __init__(self, genome_aln, insertion_aln, genome_rname, range, is_rname, clip_query_name=False):
        self.genome_aln = genome_aln
        self.insertion_aln = insertion_aln
        self.genome_rname = genome_rname
        self.range = range
        self.is_rname = is_rname

        self.clip_query_name = clip_query_name


    def run(self):
        genome_sam = pysam.AlignmentFile(self.genome_aln, "rb")
        insertion_sam = pysam.AlignmentFile(self.insertion_aln, "rb")

        genome_reads = self.get_reads_in_range(genome_sam, self.genome_rname, self.range[0], self.range[1])
        insertion_reads = self.get_insertion_reads(insertion_sam)

        out_insertion_reads = []
        genome_forward_IS_forward_count = 0
        genome_forward_IS_reverse_count = 0
        genome_reverse_IS_forward_count = 0
        genome_reverse_IS_reverse_count = 0

        for read in genome_reads:
            out_insertion_reads.append(insertion_reads[read])
            if not genome_reads[read].is_reverse and not insertion_reads[read].is_reverse:
                genome_forward_IS_forward_count += 1
            elif not genome_reads[read].is_reverse and insertion_reads[read].is_reverse:
                genome_forward_IS_reverse_count += 1
            elif genome_reads[read].is_reverse and not insertion_reads[read].is_reverse:
                genome_reverse_IS_forward_count += 1
            else:
                genome_reverse_IS_reverse_count += 1

        print 'GENOME FORWARD / IS REVERSE: %d' % genome_forward_IS_reverse_count
        print 'GENOME REVERSE / IS FORWARD: %d' % genome_reverse_IS_forward_count
        print 'GENOME FORWARD / IS FORWARD: %d' % genome_forward_IS_forward_count
        print 'GENOME REVERSE / IS REVERSE: %d' % genome_reverse_IS_reverse_count

        if genome_forward_IS_reverse_count + genome_reverse_IS_forward_count > genome_forward_IS_forward_count + genome_reverse_IS_reverse_count:
            print 'PREDICTED ORIENTATION: FORWARD'
        else:
            print 'PREDICTED ORIENTATION: REVERSE'

        outfile = pysam.AlignmentFile("allpaired.bam", "wb", template=insertion_sam)
        for read in out_insertion_reads:
            outfile.write(read)


    def get_reads_in_range(self, pysam_aln, rname, start, end):
        genome_reads = {}
        for read in pysam_aln:
            if read.reference_name == rname:
                if read.reference_start > start and read.reference_start < end:
                    query_name = read.query_name
                    if self.clip_query_name:
                        query_name = ':'.join(query_name.split(':')[:-1])

                    genome_reads[query_name] = read

        return genome_reads


    def get_insertion_reads(self, insertion_sam):
        insertion_reads = {}
        for read in insertion_sam:
            insertion_reads[read.query_name] = read
        return insertion_reads

    def count_reads(self, pysam_aln):
        count = 0
        for read in pysam_aln:
            count += 1
        return count