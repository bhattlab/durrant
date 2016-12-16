import argparse
import os, sys
from glob import glob

ref_basenames = set()

def taxonDB(path):
    files = [os.path.basename(filepath) for filepath in glob(path+'/*')]
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError('Please give absolute path to valid NCBI Taxonomy Database. Try running mustache config.')
    elif len(files) == 0:
        raise argparse.ArgumentTypeError('Please give absolute path to valid NCBI Taxonomy Database. Try running mustache config.')
    elif 'names.dmp' not in files or 'merged.dmp' not in files or 'nodes.dmp' not in files:
        raise argparse.ArgumentTypeError(
            'Please give absolute path to valid NCBI Taxonomy Database. Try running mustache config.')

    return path
