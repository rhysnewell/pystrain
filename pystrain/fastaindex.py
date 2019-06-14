import pystrain.ival as ival


class faiEntry():
    '''

    NAME	Name of this reference sequence
    LENGTH	Total length of this reference sequence, in bases
    OFFSET	Offset in the FASTA/FASTQ file of this sequence's first base
    LINEBASES	The number of bases on each line
    LINEWIDTH	The number of bytes in each line, including the newline
    QUALOFFSET	Offset of sequence's first quality within the FASTQ file
    '''

    def __init__(self, name, length, offset, linebases, linewidth, qualoffset=None):
        self.name = name
        self.length = length
        self.offset = offset
        self.linebases = linebases
        self.linewidth = linewidth
        self.qualoffset = qualoffset



class faiFile:
    """

    """

    def __init__(self, entries):
        """
        Create a vcfFile instance.
        :param entries: an iterable of entries or a filename
        """
        if isinstance(entries, str):  # filename
            self.contigs = readfaiFile(entries)
        else:
            self.contigs = dict()
            for entry in entries:
                # check if the chromosome has been seen before
                self.contigs[entry.contig] = entry

    def __len__(self):
        n = 0
        for c in self.contigs:
            n += len(self.contigs[c])
        return n

    def generate(self, contig):
        entries = self.contigs.get(contig)
        if entries != None:
            for entry in entries:
                yield entry

    def __iter__(self):
        self.contigqueue = ival.Stack()
        for contig in self.contigs.keys():
            self.contigqueue.push(self.generate(contig))
        self.current = self.contigqueue.pop()
        return self

    def __next__(self):
        try:
            ret = next(self.current)
        except StopIteration:
            if not self.contigqueue.isEmpty():
                self.current = self.contigqueue.pop()
                ret = next(self.current)
            else:
                raise StopIteration
        return ret



def readfaiFile(filename):
    """ Read a fai file
    filename: name of file
    filter_feature: name of feature to be selected, all others ignored; None means anything
    """
    contigs = dict()
    with open(filename) as f:
        row = 0
        acceptHeaderRows = 1
        headerRow = None
        for line in f:
            row += 1
            words = line.strip().split('\t')
            if len(words) == 0:
                continue  # ignore empty lines
            if words[0].strip().startswith('#'):
                continue  # comment
            if words[0].strip().startswith('browser'):
                continue  # ignore
            if words[0].strip().startswith('track'):
                continue  # ignore
            try:
                name = words[0]
                length = int(words[1])
                offset = int(words[2])

                linebases = int(words[3])
                linewidth = int(words[4])
                try:
                    qualoffset = int(words[5])
                except IndexError:
                    qualoffset = None


                entry = faiEntry(name, length, offset, linebases, linewidth, qualoffset)

                contigs[name] = entry

            except RuntimeError as e:
                if not acceptHeaderRows:
                    raise RuntimeError('Error in VCF file at row %d (%s)' % (row, e))
                else:
                    headerRow = words
                    acceptHeaderRows -= 1  # count down the number of header rows that can occur

    return contigs


if __name__ == '__main__':
    bf = faiFile('tests/r2.parent.d182.assembly.fna.fai')
    print(bf.contigs.keys())
    g = bf.generate('k141_1')
    print(next(g))
    cnt = 0
    for contig in bf:
        print(contig)