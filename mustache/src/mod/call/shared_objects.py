args = None

class pairedEndRead:

    def __init__(self, read1_in, read2_in):
        assert type(read1_in) is list, "read 1 must be a list"
        assert type(read2_in) is list, "read 2 must be a list"

        self.read1 = read1_in
        self.read2 = read2_in
        self.read2[-1] = read2_in[-1]+"\n"

        self.read1_title = read1_in[0].strip('@')
        self.read2_title = read2_in[0].strip('@')

        if len(read1_in[0]) == 0:
            raise ValueError("The input information does not represent a valid paired-"
                             "end read")

    def getReads(self):
        return [self.read1,self.read2]

    def getTitles(self):
        return [self.read1_title, self.read2_title]


class pairedEndClassification:

    def __init__(self, read1_in, read2_in):
        assert type(read1_in) is list, "read 1 must be a list"
        assert type(read2_in) is list, "read 2 must be a list"

        self.read1_title = read1_in[0].strip('@')
        self.read2_title = read2_in[0].strip('@')

        self.read1_class = read1_in[1]
        self.read2_class = read2_in[1]

        if len(read1_in[0]) == 0:
            raise ValueError("The input information does not represent a valid paired-"
                             "end read")

    def getTitles(self):
        return [self.read1_title, self.read2_title]

    def getClassifs(self):
        return [self.read1_class, self.read2_class]

def set_args(args_in):
    global args
    args = args_in