import os, subprocess
from os.path import join
version = '1.3.1'

def merge_bam_files(state, flags=['-r','-c','-p']):
    find_samtools()
    check_version()

    merged_bam_paths = {}
    for key in state.paths.bam_info:
        bamfiles = [state.paths.bam_info[key][sample] for sample in state.paths.bam_info[key]]
        outbam = join(state.paths.merged_bam_dir, key+'.bam')
        command = ['samtools','merge'] + flags + [join(state.paths.merged_bam_dir, key+'.bam')] + bamfiles
        merged_bam_paths[key] = outbam
        subprocess.check_output(command)

    state.paths.merged_bam_paths = merged_bam_paths

def process_all_taxon_sams(state):
    find_samtools()
    check_version()

    merged_bam_paths = {}
    for key in state.paths.merged_bam_paths:
        path = state.paths.merged_bam_paths[key]
        bamfile = bamsort(path)
        bamfile = bamindex(bamfile)
        merged_bam_paths[key] = bamfile

    state.paths.merged_bam_paths = merged_bam_paths

def sam_to_bam(samfile):
    prefix = os.path.join(os.path.dirname(samfile), '.'.join(os.path.basename(samfile).split('.')[:-1]))
    bam_out = prefix+'.bam'
    samtools_args = " ".join(['samtools view -b -S -o', bam_out, samfile])

    subprocess.call(samtools_args, shell=True)
    return bam_out

def bamsort(bamfile):
    prefix = os.path.join(os.path.dirname(bamfile), os.path.basename(bamfile).split('.')[0])
    bam_out = prefix+".bam"
    samtools_args = "samtools sort -T %s -o %s %s" % (prefix+".sorted", bam_out, prefix+".bam")
    subprocess.call(samtools_args, shell=True)
    return bam_out

def bamindex(bamfile):
    samtools_args = "samtools index %s" % bamfile

    subprocess.call(samtools_args, shell=True)
    return bamfile

def find_samtools():
    try:
        subprocess.check_output(['samtools', '--help'])
    except OSError:
        raise RuntimeError('samtools not found; check if it is installed and in $PATH\n')

def check_version():
    global version

    stdout = subprocess.check_output(['samtools', '--version'])
    local_version = stdout.split('\n')[0].split()[-1]
    assert version == local_version, 'samtools version incompatibility %s != %s' % (version, local_version)