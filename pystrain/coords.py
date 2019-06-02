#!/usr/bin/env python3
import string, sys, re, math, os, array
import numpy
import pystrain.ival as ival
import pystrain.sequence as sequence


class nucmerCoords(object):
    s1_start = None
    s1_end = None
    s1_len = None
    s1_strand = None
    s2_start = None
    s2_end = None
    s2_len = None
    s2_strand = None
    percent_id = None
    r_len = None
    r_cov = None
    r_tag = None
    q_len = None
    q_cov = None
    q_tag = None

    def __init__(self, entry):
        self.s1_start = int(entry[0])
        self.s1_end = int(entry[1])
        if self.s1_start > self.s1_end:
            self.s1_strand = '-'
            self.s1_start = int(entry[1])
            self.s1_end = int(entry[0])
        else:
            self.s1_strand = '+'
        self.s1_len = int(entry[4])
        self.s2_start = int(entry[2])
        self.s2_end = int(entry[3])
        if self.s2_start > self.s2_end:
            self.s2_strand = '-'
            self.s2_start = int(entry[3])
            self.s2_end = int(entry[2])
        else:
            self.s2_strand = '+'
        self.s2_len = int(entry[5])
        self.percent_id = float(entry[6])
        self.r_len = int(entry[7])
        self.r_cov = float(entry[9])
        self.r_tag = str(entry[11])
        self.q_len = int(entry[8])
        self.q_cov = float(entry[10])
        self.q_tag = str(entry[12])

    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        entry = format("%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%.3f\t%s\t%s" %
                       (self.s1_start, self.s1_end, self.s1_strand, self.s2_start, self.s2_end, self.s2_strand,
                        self.s1_len, self.s2_len, self.percent_id, self.r_len, self.q_len, self.r_cov, self.q_cov,
                        self.r_tag, self.q_tag))
        return entry


class coordFile():

    def __init__(self, entries, reference_path, query_path):
        self.source = dict()
        self.contigs = dict()
        self.reference = sequence.readFastaFile(reference_path)
        self.query = sequence.readFastaFile(query_path)
        for entry in entries:
            # check if the ref contig has been seen before
            ref_tree = self.contigs.get(entry.r_tag)
            if not ref_tree:
                ref_tree = ival.IntervalTree()
                self.contigs[entry.r_tag] = ref_tree
                self.source[entry.r_tag] = 'r'
            # put the entry in the reference interval tree for the appropriate contig
            iv = ival.Interval(entry.s1_start, entry.s1_end)
            ref_tree.put(iv, entry)

            # check if the query contig has been seen before
            query_tree = self.contigs.get(entry.q_tag)
            if not query_tree:
                query_tree = ival.IntervalTree()
                self.contigs[entry.q_tag] = query_tree
                self.source[entry.q_tag] = 'q'
            # put the entry in the interval tree for the appropriate contig
            iv = ival.Interval(entry.s2_start, entry.s2_end)
            query_tree.put(iv, entry)

    def __len__(self):
        n = 0
        for c in self.contigs:
            n += len(self.contigs[c])
        return n

    def generate(self, contig):
        mytree = self.contigs.get(contig)
        if mytree is not None:
            for e in mytree:
                for entry in e.values:
                    yield entry

    def __iter__(self):
        self.contigqueue = ival.Stack()
        for c in sorted(self.contigs.keys())[::-1]:
            self.contigqueue.push(self.generate(c))
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
        if isinstance(item, nucmerCoords):
            tree = self.contigs.get(item.r_tag)
            if tree == None: return False
            else: return ival.Interval(item.s1_start, item.s1_end) in tree
        else:
            return False

    def getOverlap(self, item):
        """
        :param item: a nucmerCoord entry
        :return: an array of lists of overlaps, first item is reference overlaps, second is query overlaps
        """
        if isinstance(item, nucmerCoords):
            overlaps = {}
            tree = self.contigs.get(item.r_tag)
            if tree == None:
                overlaps[item.r_tag] = [None]
            else:
                iv = ival.Interval(item.s1_start, item.s1_end)
                res = tree.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                overlaps[item.r_tag] = ret

            tree = self.contigs.get(item.q_tag)
            if tree == None:
                overlaps[item.q_tag] = [None]
            else:
                iv = ival.Interval(item.s2_start, item.s2_end)
                res = tree.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                overlaps[item.q_tag] = ret

            return overlaps
        else:
            return {'reference':[None], 'query':[None]}

    def getClosest(self, item):
        if isinstance(item, nucmerCoords):
            tree = self.contigs.get(item.r_tag)
            if tree == None: return None
            else:
                iv = ival.Interval(item.s1_start, item.s1_end)
                node = tree.closest(iv)
                if node != None: return node.values
                else: return None
        else: return None

    def getOneOfClosest(self, item):
        all = self.getClosest(item)
        if all == None: return None
        else: return next(iter(all))

    def getOneOfOverlap(self, item):
        all = self.getOverlap(item)
        if all == None: return None
        elif len(all) == 0: return None
        else: return next(iter(all))

def readCoordFile(filename):

    fh = open(filename)
    coordlist = []
    batch = '' # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    start = False
    file_line = True
    ref_path = None
    query_path = None
    for row in fh:
        row = row.strip()
        if file_line:
            row = row.strip().split()
            ref_path = row[0]
            query_path = row[1]
            file_line = False

        elif row.startswith('===='):
            start = True
        elif start is True:
            row = ''.join(row.split('|')).split()
            coordlist.append(nucmerCoords(row))
    fh.close()
    coords_file = coordFile(coordlist, ref_path, query_path)
    return coords_file

if __name__=='__main__':
    filename = sys.argv[1]
    print(filename)
    coords = readCoordFile(filename)
    print(coords)
    for coord in coords:
        print(coord)