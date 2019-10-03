import pystrain.ival as ival
from collections import OrderedDict
import numpy as np
from sklearn.preprocessing import scale, MinMaxScaler

class contigEntry():

    def __init__(self, fields):
        self.contigName = fields[0]
        self.contigLen = int(fields[1])
        self.totalAvgDepth = float(fields[2])
        self.totalAvgGeno = float(fields[3])
        self.sampleDepths = []
        self.sampleVars = []
        self.sampleGenos = []
        sampleFields = [fields[pos:pos+3] for pos in range(4, len(fields), 3)]
        for field in sampleFields:
            self.sampleDepths.append(float(field[0]))
            self.sampleVars.append(float(field[1]))
            self.sampleGenos.append(float(field[2]))

    def __str__(self):
        entry = "{}\t{}\t{}\t{}".format(self.contigName, self.contigLen, self.totalAvgDepth, self.totalAvgGeno)
        values = ""
        for c, v, g in zip(self.sampleDepths, self.sampleVars, self.sampleGenos):
            values += "\t{}\t{}\t{}".format(c, v, g)

        entry = entry + values
        return entry

    def extract(self):
        entry = [self.contigLen, self.totalAvgDepth, self.totalAvgGeno]
        for c, v, g in zip(self.sampleDepths, self.sampleVars, self.sampleGenos):
            entry.append(c)
            entry.append(v)
            entry.append(g)
        return entry




class contigStats():

    def __init__(self, entries):
        """
        Create a contigStats instance.
        :param entries: an iterable of entries or a filename
        """
        if isinstance(entries, str):  # filename
            self.contigs = readstatsFile(entries)
        else:
            self.contigs = OrderedDict()
            for entry in entries:
                self.contigs[entry.contigName] = entry

    def __len__(self):
        return len(self.contigs)

    def generate(self, contig):
        entries = self.contigs.get(contig)
        return entries

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

    def array(self, transpose=True, minmax=True):
        # colvals = np.array(list(self.contigs.values()))
        colvals = [0]*len(self.contigs)
        for idx, entry in enumerate(self.contigs.values()):
            colvals[idx] = entry.extract()
        colvals = np.array(colvals)
        if minmax is True:
            scaler = MinMaxScaler(feature_range=(1,2))
            scaler.fit(colvals)
            colvals = scaler.transform(colvals)
            colvals = np.nan_to_num(colvals)

        else:
            colvals = scale(colvals, axis=0)
            colvals = np.nan_to_num(colvals)

        if transpose is True:
            colvals = colvals.T


        return colvals





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