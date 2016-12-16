import sys
from collections import defaultdict, OrderedDict


def read_genome_alignment(sam_file, start_read):

    sam_file, first_line = read_genome_alignment_fast(sam_file, start_read)
    sam_header = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
    current_set = None
    while True:
        try:
            first_line[-1] = first_line[-1].strip()
            first_line = OrderedDict(zip(sam_header, first_line))
            current_name = first_line['QNAME']
            current_set = [first_line]

            for line in sam_file:
                line = OrderedDict(zip(sam_header, line.strip().split()))

                name = line['QNAME']
                if name == current_name:
                    current_set.append(line)
                else:
                    yield current_set
                    current_name = name
                    current_set = [line]

        except KeyError:
            break

    yield current_set

def read_genome_alignment_fast(sam_file, start_read):
    read_count= 1
    while True:

        first_line = sam_file.readline()
        while first_line[0] == '@': first_line = sam_file.readline()
        first_line = first_line.split()
        if start_read == 1: return sam_file, first_line

        current_name = first_line[0]
        current_set = [first_line]

        for line in sam_file:
            line = line.split()
            name = line[0]
            if name == current_name:
                current_set.append(line)
            else:
                current_name = name
                current_set = [line]
                read_count += 1
                if read_count == start_read:
                    return sam_file, line

    raise IndexError

def read_insertion_alignments(sam_file, max_lines):
    sam_file, first_line = read_insertion_alignments_fast(sam_file, max_lines)
    sam_header = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
    current_set = None

    while True:
        try:
            first_line[-1] = first_line[-1].strip()
            first_line = OrderedDict(zip(sam_header, first_line))
            current_name = first_line['QNAME']
            current_set = [first_line]

            for line in sam_file:
                line = dict(zip(sam_header, line.strip().split()))

                name = line['QNAME']
                if name == current_name:
                    current_set.append(line)
                else:
                    yield current_set
                    current_name = name
                    current_set = [line]

        except KeyError:
            break

    yield current_set

def read_insertion_alignments_fast(sam_file, start_read):
    read_count= 1
    while True:

        first_line = sam_file.readline()
        while first_line[0] == '@': first_line = sam_file.readline()
        first_line = first_line.split()
        if start_read == 1: return sam_file, first_line

        current_name = first_line[0]
        current_set = [first_line]

        for line in sam_file:
            line = line.split()
            name = line[0]
            if name == current_name:
                current_set.append(line)
            else:
                current_name = name
                current_set = [line]
                read_count += 1
                if read_count == start_read:
                    return sam_file, line

    raise IndexError("Start read is greater than file size")

def read_classifications(read_class_file, max_lines, delim='\t'):
    read_class_file, header = read_classifications_fast(read_class_file, max_lines)
    while True:
        try:
            for line in read_class_file:
                read_class = dict(zip(header, line.strip().split('\t')))
                read_class['HEADER'] = read_class['HEADER'].split()[0].strip('@')
                yield read_class
        except ValueError:
            break

def read_classifications_fast(read_class_file, start_read, delim='\t'):
    header = ''
    header = read_class_file.readline().strip().split(delim)
    for title in range(len(header)):
        header[title] = header[title].replace(' ', '').upper()

    line_count = 0
    if line_count == start_read-1:
        return read_class_file, header

    for line in read_class_file:
        line_count += 1
        if line_count == start_read-1:
            return read_class_file, header

    raise IndexError("Start read is greater than file size")