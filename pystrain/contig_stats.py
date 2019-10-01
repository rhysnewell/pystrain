import pystrain.ival as ival


class contigEntry():

    def __init__(self, fields):
        self.contigName = fields[0]
        self.contigLen = fields[1]
        self.totalAvgDepth = fields[2]
        self.totalAvgGeno = fields[3]
        self.sampleDepths = []
        self.sampleVars = []
        self.sampleGenos = []
        sampleFields = [fields[pos:pos+3] for pos in range(4, len(fields), 3)]
        for field in sampleFields:
            self.sampleDepths.append(field[0])
            self.sampleVars.append(field[1])
            self.sampleGenos.append(field[2])


class contigStats():

    def __init__(self, entries):
        """
        Create a contigStats instance.
        :param entries: an iterable of entries or a filename
        """
        if isinstance(entries, str):  # filename
            self.contigs = readstatsFile(entries)
        else:
            self.contigs = dict()
            for entry in entries:
                self.contigs[entry.contigName] = entry

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

    def __contains__(self, item):
        if isinstance(item, contigEntry):
            tree = self.contigs.get(item.contigName)
            if tree == None: return False
            else: return True
        else:
            return False




def readstatsFile(filename):
    """ Read a contig stats file from lorikeet binning output.
    filename: name of file
    filter_feature: name of feature to be selected, all others ignored; None means anything
    """
    contigs = dict()
    with open(filename) as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            else:
                line = line.strip().split()
                entry = contigEntry(line)
                contigs[entry.contigName] = entry
    return contigs