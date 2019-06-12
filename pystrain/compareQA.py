#!/usr/bin/env python3

import sys, os
class qaResult(object):
    bin = None
    completeness = None
    contamination = None
    heterogeneity = None


    def __init__(self, entry):
        self.bin = entry[0].strip('new_')
        self.completeness = float(entry[11])
        self.contamination = float(entry[12])
        self.heterogeneity = float(entry[13])


    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        entry = format("%s\t%.3f\t%.3f\t%.3f" %
                       (self.bin, self.completeness, self.contamination, self.heterogeneity))
        return entry


class qaFile():

    def __init__(self, entries):
        self.bins = dict()

        for entry in entries:
            self.bins[entry.bin] = entry


    def __len__(self):
        n = len(self.bins.keys())
        return n




def readCheckMQAFile(filename):
    qalist = []
    start = 0

    with open(filename) as fh:
        for row in fh:
            row = row.strip()
            if start == 2:
                row = ''.join(row.split('|')).split()
                qalist.append(qaResult(row))
            elif row.startswith('------'):
                start += 1

    qa_file = qaFile(qalist)
    return qa_file


def compareQA(original, new):
    print("bin\tcompleteness_increase\tcontamination_increase\theterogenity_increase")
    for bin in original.bins.keys():
        if bin in new.bins.keys():
            print("%s\t%.3f\t%.3f\t%.3f" % (bin, (new[bin].completeness - original[bin].completeness),
                                            (new[bin].contamination - original[bin].contamination),
                                            (new[bin].heterogeneity - original[bin].heterogeneity)))

if __name__=="__main__":
    try:
        original_path = sys.argv[1]
        new_path = sys.argv[2]
        original = readCheckMQAFile(original_path)
        new = readCheckMQAFile(new_path)
        compareQA(original, new)
    except IndexError:
        print("Usage: <Original checkM QA> <New checkM QA> \n"
              "compareQA.py used to compare checkM results on sets of bins. Only works if bin names correspond"
              "between assemblies")