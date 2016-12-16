import IO
import peaks


def test_readers(gen1, gen2, class1, class2, is1, is2):
    line_start = 23
    gen1_in = IO.read_genome_alignment(open(gen1), line_start)
    gen2_in = IO.read_genome_alignment(open(gen2), line_start)

    is1_in = IO.read_insertion_alignments(open(is1), line_start)
    is2_in = IO.read_insertion_alignments(open(is2), line_start)

    class1_in = IO.read_classifications(open(class1), line_start)
    class2_in = IO.read_classifications(open(class2), line_start)
    print gen1_in.next()
    print gen2_in.next()

    print is1_in.next()
    print is2_in.next()

    print class1_in.next()
    print class2_in.next()


if __name__ == "__main__":
    #gen1 = "data/fastq1_to_genome_bcaccae_sequence.sam"
    #gen2 = "data/fastq2_to_genome_bcaccae_sequence.sam"
    #is1 = "data/fastq1_to_IS.sam"
    #is2 = "data/fastq2_to_IS.sam"
    #class1 = "data/small_PE1.fq.results.tsv"
    #class2 = "data/small_PE1.fq.results.tsv"

    #test_readers(gen1, gen2, class1, class2, is1, is2)

    #genome1 = "/Users/mdurrant/OneDrive/Stanford/BhattLab/code/bhattlab/mytools/ISPeaks/data/refs/bcaccae_sequence.fasta"
    #genome2 = "/Users/mdurrant/OneDrive/Stanford/BhattLab/code/bhattlab/mytools/ISPeaks/data/IS_fastas/Bacteroides_all.fasta"
    #print misc.get_genome_length(genome1)
    #print misc.get_genome_length(genome2)

    peaks1 = [[1,100],[50,150]]
    peaks2 = [[1,100],[101,200]]
    peaks3 = [[1, 100], [101, 200], [201, 300], [301, 100000]]
    peaks4 = [[1, 100], [99, 100], [150,152], [301, 100000],[500,100001]]
    peaks5 = [[2,100],[3,1000],[4,500], [1,100000]]


    print "IN:", peaks2
    print "OUT:", peaks.merge_peaks_single(peaks2)

    print
    print "IN:", peaks2
    print "OUT:", peaks.merge_peaks_single(peaks2, 1)

    print
    print "IN:", []
    print "OUT:", peaks.merge_peaks_single([], 1)

    print
    print "IN:", [[1,100]]
    print "OUT:", peaks.merge_peaks_single([[1, 100]], 1)

    print
    print "IN:", [[1, 100]]
    print "OUT:", peaks.merge_peaks_single([[1, 100]], 1)

    print
    print "IN:", peaks3
    print "OUT:", peaks.merge_peaks_single(peaks3, 0)

    print
    print "IN:", peaks3
    print "OUT:", peaks.merge_peaks_single(peaks3, 1)

    print
    print "IN:", peaks3
    print "OUT:", peaks.merge_peaks_single(peaks3, 1000)

    print
    print "IN:", peaks4
    print "OUT:", peaks.merge_peaks_single(peaks4)

    print
    print "IN:", peaks5
    print "OUT:", peaks.merge_peaks_single(peaks5)

    print
    print "IN:", peaks1, peaks2
    print "OUT:", peaks.merge_peaks_double(peaks1, peaks2)