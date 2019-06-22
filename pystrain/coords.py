#!/usr/bin/env python3
import string, sys, re, math, os, array
import numpy
import pyfaidx
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
                        self.r_tag.strip(), self.q_tag.strip()))
        return entry


class coordFile():

    def __init__(self, entries, reference_path=None, query_path=None):
        # self.source = dict()
        # self.contigs = dict()
        self.r_contigs = dict()
        self.q_contigs = dict()
        if reference_path is not None:
            self.reference_name = reference_path.split('/')[-1]
            self.reference = pyfaidx.Faidx(reference_path)
        if query_path is not None:
            self.query_name = query_path.split('/')[-1]
            self.query = pyfaidx.Faidx(query_path)


        for entry in entries:
            # check if the ref contig has been seen before
            ref_tree = self.r_contigs.get(entry.r_tag)
            if not ref_tree:
                ref_tree = ival.IntervalTree()
                self.r_contigs[entry.r_tag] = ref_tree
                # self.source[entry.r_tag] = 'r'
            # put the entry in the reference interval tree for the appropriate contig
            iv = ival.Interval(entry.s1_start, entry.s1_end)
            ref_tree.put(iv, entry)

            # check if the query contig has been seen before
            query_tree = self.q_contigs.get(entry.q_tag)
            if not query_tree:
                query_tree = ival.IntervalTree()
                self.q_contigs[entry.q_tag] = query_tree
                # self.source[entry.q_tag] = 'q'
            # put the entry in the interval tree for the appropriate contig
            iv = ival.Interval(entry.s2_start, entry.s2_end)
            query_tree.put(iv, entry)

    def __len__(self):
        n = 0
        for c in self.r_contigs:
            n += len(self.r_contigs[c])
        return n

    def generate(self, contig, source='r'):
        """

        :param contig: The contig name
        :param source: Whether you want to reference or query contig. Contig names can be shared across assemblies
        :return: A coordinate for that contig
        """
        if source.lower().startswith('r'):
            mytree = self.r_contigs.get(contig)
            if mytree is not None:
                for e in mytree:
                    for entry in e.values:
                        yield entry
        elif source.lower().startswith('q'):
            mytree = self.q_contigs.get(contig)
            if mytree is not None:
                for e in mytree:
                    for entry in e.values:
                        yield entry

    def __iter__(self):
        self.contigqueue = ival.Stack()
        for c in sorted(self.r_contigs.keys())[::-1]:
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
            tree = self.r_contigs.get(item.r_tag)
            if tree == None: return False
            else: return ival.Interval(item.s1_start, item.s1_end) in tree
        else:
            return False


def readCoordFile(filename, readFastas=True):
    coordlist = []
    batch = ''  # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    start = False
    file_line = True
    ref_path = None
    query_path = None
    with open(filename) as fh:
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
    if readFastas:
        coords_file = coordFile(coordlist, ref_path, query_path)
    else:
        coords_file = coordFile(coordlist)
    return coords_file

if __name__=='__main__':
    filename = sys.argv[1]
    print(filename)
    coords = readCoordFile(filename)
    print(coords)
    for coord in coords.generate('query_contig1', source='q'):
        print(coord)