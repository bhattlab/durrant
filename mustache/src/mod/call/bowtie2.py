"""
Wrapper script for bowtie2-align and bowtie2-build with version control
"""
import subprocess
from glob import glob
import os, sys

args = None

def build(version, fasta, logger=None):
    global args
    if not logger:
        if args: logger = args['logger']
    """
    Construct bowtie2 indices (.bt2 files)
    :param version: Enforces bowtie2 version number (e.g., '2.2.1')
    :param fasta: Path to FASTA containing reference sequences
    :return:
    """
    find_bowtie2()
    check_version(version)
    if is_built(fasta):
        if logger: logger.info("The fasta appears to have already been built and indexed...")
        return

    stdout = subprocess.check_output(['bowtie2-build', '--version'])
    local_version = stdout.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie2-build version incompatibility %s != %s' % (version, local_version)
    subprocess.check_call(['bowtie2-build', '--offrate', '4', '-q', '-f', fasta, fasta])


def align_fastq_to_insertions(version, refpath, fastq, outfile, threads=1,
                              flags=('--local', '--quiet', '--reorder','--all')):
    # check that we have access to bowtie2
    find_bowtie2()

    # check that version is the expected version
    check_version(version)

    # stream output from bowtie2
    bowtie_args = " ".join(['bowtie2', '-x', refpath, '-U', fastq, '-S %s' %
                            outfile, '-p', str(threads)] + list(flags))

    subprocess.call(bowtie_args, shell=True)
    return outfile

def align_to_genome(version, refpath, fastq, outfile, threads=1,
                              flags=('--local', '--quiet', '--reorder','--all')):
    # check that we have access to bowtie2
    find_bowtie2()

    # check that version is the expected version
    check_version(version)

    # stream output from bowtie2
    bowtie_args = " ".join(['bowtie2', '-x', refpath, '-U', fastq, '-S %s' %
                             outfile, '-p', str(threads)] + list(flags))

    subprocess.call(bowtie_args, shell=True)
    return outfile

def find_bowtie2():
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

def check_version(version):
    stdout = subprocess.check_output(['bowtie2', '--version'])
    local_version = stdout.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie2 version incompatibility %s != %s' % (version, local_version)

def is_built(fasta):
    if len(glob(fasta+"*")) > 1:
        return True
    return False

def sam_to_bam(samfile):
    prefix = os.path.join(os.path.dirname(samfile), os.path.basename(samfile).split('.')[0])
    bam_out = prefix+'.bam'
    samtools_args = " ".join(['samtools view -b -S -o', bam_out, samfile])

    subprocess.call(samtools_args, shell=True)
    return bam_out

def bamsort(bamfile):
    prefix = os.path.join(os.path.dirname(bamfile), os.path.basename(bamfile).split('.')[0])
    bam_out = prefix+".sorted.bam"
    samtools_args = "samtools sort -T %s -o %s %s" % (prefix+".sorted", bam_out, prefix+".bam")
    subprocess.call(samtools_args, shell=True)
    return bam_out

def bamindex(bamfile):
    samtools_args = "samtools index %s" % bamfile

    subprocess.call(samtools_args, shell=True)
    return bamfile