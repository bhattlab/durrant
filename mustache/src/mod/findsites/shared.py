

def get_peaks(peaks_path):
    peaks_out = []
    with open(peaks_path) as peaks:
        for peak in peaks:
            peak = peak.strip().split('\t')
            peaks_out.append([peak[0], int(peak[1]) - 1, int(peak[2]) - 1])

    return peaks_out


## Miscellaneous
def peak_str(start, end, in_zeroi=True, out_zeroi=False):

    if in_zeroi:
        if out_zeroi:
            return '-'.join([str(start), str(end)])
        else:
            return '-'.join([str(start + 1), str(end + 1)])
    else:
        if out_zeroi:
            return '-'.join([str(start - 1), str(end - 1)])
        else:
            return '-'.join([str(start), str(end)])

def pprint_list(list_in):
    return ', '.join([str(elem) for elem in list_in])