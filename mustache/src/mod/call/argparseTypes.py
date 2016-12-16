import argparse
import os

import misc

ref_basenames = set()

def genome(path):
    global ref_basenames
    try:
        if not os.path.isfile(path): raise TypeError("Genome fasta is not a file.")
        if misc.get_refid_from_path(path) in ref_basenames: raise TypeError("Duplicate Basename")
        ref_basenames.add(os.path.basename(path).split('.')[0])
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError("Genomes must be entered as \"<reference-fasta1> <reference-fasta2>\" and "
                                         "they must have unique basenames.")

def insertion_fasta(path):
    try:
        if not os.path.isfile(path): raise TypeError()
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Please give a path to a valid fasta file.')


def fastq_file(path):
    try:
        if not os.path.isfile(path): raise TypeError()
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Please give paths to valid fastq files.')


def class_file(path):
    try:
        if not os.path.isfile(path): raise TypeError()

        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Please give paths to valid classification files. Make sure that it contains'
                                         ' an alphanumeric characters and white space delimiters only.')

def complete_class_exclude(exclude):
    try:
        if '=' not in exclude: raise TypeError()
        if len(exclude.split('=')) != 2: raise TypeError()

        exclude = exclude.split('=')[0].upper(), exclude.split('=')[1]
        return exclude
    except:
        raise argparse.ArgumentTypeError('--classification-exclude must be formatted as \"<Column Name1>=<Criterion1> <Column Name2>=<Criterion2>...\"')

def genome_class_exclude(exclude):
    try:
        if '=' not in exclude: raise TypeError()
        if len(exclude.split('=')) != 2: raise TypeError()

        exclude = exclude.split('=')[0].upper(), exclude.split('=')[1]
        return exclude
    except:
        raise argparse.ArgumentTypeError('--classification-exclude must be formatted as \"<Column Name1>=<Criterion1> <Column Name2>=<Criterion2>...\"')

def IS_class_exclude(exclude):
    try:
        if '=' not in exclude: raise TypeError()
        if len(exclude.split('=')) != 2: raise TypeError()

        exclude = exclude.split('=')[0].upper(), exclude.split('=')[1]
        return exclude
    except:
        raise argparse.ArgumentTypeError(
            '--classification-exclude must be formatted as \"<Column Name1>=<Criterion1> <Column Name2>=<Criterion2>...\"')

def output_folder(path):
    try:
        if os.path.isdir(path):
            raise TypeError()
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Output folder cannot already exist')


def taxon_nodes(path):
    try:
        if not os.path.isfile(path): raise TypeError()
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Please give paths to valid taxonomy node files.')

def taxon_names(path):
    try:
        if not os.path.isfile(path): raise TypeError()
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Please give paths to valid taxonomy names file.')
