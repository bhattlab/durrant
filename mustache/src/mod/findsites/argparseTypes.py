import argparse
import os, sys
import pysam

ref_basenames = set()

def bam_file(path):

    try:
        if not os.path.isfile(path): raise TypeError()

        handle = pysam.AlignmentFile(path, 'rb')
        read = handle.next()
        return os.path.abspath(path)

    except:
        raise argparse.ArgumentTypeError('Please gave path of valid bam file.')

def peak_ranges(path):

    try:
        if not os.path.isfile(path): raise TypeError()
        with open(path, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                if len(line) != 3: raise TypeError()
                start = int(line[1])
                stop = int(line[2])
                if start >= stop: raise TypeError()

        return os.path.abspath(path)

    except:
        raise argparse.ArgumentTypeError('Peak Ranges file is not in the correct format.')


def genome(path):
    global ref_basenames
    try:
        if not os.path.isfile(path): raise TypeError("Genome fasta is not a file.")
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError("Genomes must be entered fasta file.")
