import subprocess
import sys
from glob import glob
import misc, samtools
import pysam
from os.path import basename, join

def call_peaks(version, bampath, genome_length, outdir):
    # check that we have access to macs2
    find_macs2()

    # check that version is the expected version
    # check_version(version)
    # stream output from bowtie2
    prefix_out = ".".join(basename(bampath).split('.')[:-1])
    macs2_args = 'macs2 callpeak -t %s -g %d -n %s --broad --broad-cutoff 0.1 --nomodel --extsize 147 --outdir %s --verbose 0' \
                 % (bampath, genome_length, prefix_out, outdir)

    subprocess.call(macs2_args.split(), shell=False)
    broadPeak_file = glob(join(outdir, prefix_out)+'_*broadPeak')

    if len(broadPeak_file) != 1:
        raise RuntimeError("Figure this out. Be a better programmer.")
    return broadPeak_file[0]


def find_macs2():
    try:
        subprocess.check_output(['macs2', '-h'])
    except OSError:
        raise RuntimeError('macs2 not found; check if it is installed and in $PATH\n')