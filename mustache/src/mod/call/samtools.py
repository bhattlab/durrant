import os, subprocess

version = '1.3.1'

def process_all_taxon_sams(state):
    find_samtools()
    check_version()

    for key in state.paths.taxon_sam_paths:
        path = state.paths.taxon_sam_paths[key]
        bamfile = sam_to_bam(path)
        bamfile = bamsort(bamfile)
        bamfile = bamindex(bamfile)
        state.paths.taxon_bam_paths[key] = bamfile

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