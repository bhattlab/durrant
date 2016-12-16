import argparse
import os, sys

ref_basenames = set()

def ispeaks_directory(path):

    try:
        if not os.path.isdir(path): raise TypeError()
        return os.path.abspath(path)
    except:
        raise argparse.ArgumentTypeError('Please give paths to valid ispeaks output directories files.')

def output_folder(path):
    try:
        if os.path.isdir(path):
            raise TypeError()
        return os.path.abspath(path)

    except:
        raise argparse.ArgumentTypeError('Output folder cannot already exist')


def minimum_length(nmin):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values):
                msg='argument "{f}" requires at least {nmin} arguments'.format(
                    f=self.dest,nmin=nmin)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength