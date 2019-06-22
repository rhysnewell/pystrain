"""
Module *** sequence ***

This module depends on the following modules

sym -- defines an alphabet
prob -- defines structures to hold probabilities (prob also depends on sym)

This module incorporates classes for

Sequence -- names and defines a sequence of symbols; computes various transformations and pairwise alignments
Alignment -- defines a multiple sequence alignment; computes stats for use in substitution matrices
SubstMatrix -- substitution matrix class to support alignment methods
Regexp -- defines patterns as regular expressions for textual pattern matching in sequences
PWM -- defines a weight matrix that can score any site in actual sequences

Incorporates methods for loading and saving files relevant to the above (e.g. FASTA, ALN, substitution matrices)
and methods for retrieving relevant data from web services

This code has been adapted to Python 3.5 in 2017

This code has gone through many updates and has benefited from kind contributions of course participants.
Please keep suggestions coming!
Email: m.boden@uq.edu.au
"""

import string, sys, re, math, os, array
from pystrain.prob import *
import pystrain.ival as ival

# Sequence ------------------****

class Sequence(object):
    """ A biological sequence. Stores the sequence itself (as a compact array), 
    the alphabet (i.e., type of sequence it is), and optionally a name and further 
    information. """
    
    sequence = None # The array of symbols that make up the sequence 
    alphabet = None # The alphabet from which symbols come
    name =     None # The name (identifier) of a sequence
    info =     None # Other information (free text; e.g. annotations)
    length =   None # The number of symbols that the sequence is composed of
    gappy =    None # True if the sequence has "gaps", i.e. positions that represent deletions relative another sequence
    fragmentStart = None
    fragmentEnd = None
    percentIdentity = None
    
    def __init__(self, sequence, alphabet = None, name = '', info = '',
                 gappy = False, fragmentposition = None, percentid = None):
        """ Create a sequence with the sequence data. Specifying the alphabet,
        name and other information about the sequence are all optional.
        The sequence data is immutable (stored as a string).
        Example:
        >>> myseq = Sequence('MVSAKKVPAIAMSFGVSF')
        will create a sequence with no name, and assign one of the predefined
        alphabets on the basis of what symbols were used.
        >>> myseq.alphabet.symbols
        will output the standard protein alphabet:
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
        'R', 'S', 'T', 'V', 'W', 'Y']
        fragment specifies whether this sequence is a fragment of a larger sequence e.g. a contig and its input
        should be a list with the start and end positions of the fragment e.g. [250, 1000]"""
        
        self.sequence = sequence
        
        # Assign an alphabet
        # If no alphabet is provided, attempts to identify the alphabet from sequence
        self.alphabet = None
        if not alphabet is None:
            for sym in self.sequence:
                if not sym in alphabet and (sym != '-' or not gappy):  # error check: bail out
                    raise RuntimeError('Invalid symbol: %c in sequence %s' % (sym, name))
            self.alphabet = alphabet
        else:
            for alphaName in preferredOrder:
                alpha = predefAlphabets[alphaName]
                valid = True
                for sym in self.sequence:
                    if not sym in alpha and (sym != '-' or not gappy):  
                        valid = False
                        break
                if valid:
                    self.alphabet = alpha
                    break
            if self.alphabet is None:
                raise RuntimeError('Could not identify alphabet from sequence: %s' % name)


        # Store other information
        self.name = name
        self.info = info
        self.length = len(self.sequence)
        self.gappy = gappy

        #define which fragment of a contig we have
        if fragmentposition is None:
            self.fragmentStart = 1
            self.fragmentEnd = self.length
        else:
            self.fragmentStart = fragmentposition[0]
            self.fragmentEnd = fragmentposition[1]

        if percentid is None:
            self.percentIdentity = 100
        else:
            self.percentIdentity = percentid
        
    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print (len(seq))
        9
        """
        return len(self.sequence)

    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        str = self.name + ': '
        for sym in self:
            str += sym
        return str
    
    def __iter__(self):
        """ Defines how a Sequence should be "iterated", i.e. what its elements are, e.g.
        >>> seq = Sequence('AGGAT', DNA_Alphabet)
        >>> for sym in seq:
                print (sym)
        will print A, G, G, A, T (each on a separate row)
        """ 
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()
    
    def __contains__(self, item):
        """ Defines what is returned when the "in" operator is used on a Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print ('T' in seq)
        True
            which is equivalent to 
        >>> print (seq.__contains__('T'))
        True
        >>> print ('X' in seq)
        False
        """ 
        for sym in self.sequence:
            if sym == item:
                return True
        return False
        
    def __getitem__(self, ndx):
        """ Retrieve a specified index (or a "slice" of indices) of the sequence data.
            Calling self.__getitem__(3) is equivalent to self[3] 
        """
        if type(ndx) is slice:
            return ''.join(self.sequence[ndx])
        else:
            return self.sequence[ndx]
        
    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        if parseDefline(self.info)[0] == self.name: # this sequence was previously "parsed" and info should hold the original header
            fasta = '>' + self.info + '\n'
        else:
            fasta = '>' + self.name.strip('_rq') + ' ' + self.info + '\n'
        fasta += ''.join(self.sequence)
        fasta += '\n'

        return fasta
    
    def count(self, findme = None):
        """ Get the number of occurrences of specified symbol findme OR
            if findme = None, return a dictionary of counts of all symbols in alphabet """
        if findme != None:
            cnt = 0
            for sym in self.sequence:
                if findme == sym:
                    cnt = cnt + 1
            return cnt
        else:
            symbolCounts = {}
            for symbol in self.alphabet:
                symbolCounts[symbol] = self.count(symbol)
            return symbolCounts

    def getDegapped(self):
        """ Create the sequence excluding gaps, and provide the corresponding indices for the gapped version, e.g.
        >>> gappy = Sequence('AC--TA-GA', DNA_Alphabet, name = 'myseq', gappy = True)
        >>> degapped, indices = gappy.getDegapped()
        >>> print(degapped)
            myseq: ACTAGA
        >>> print(indices)
            [0, 1, 4, 5, 7, 8]
        """
        idxs = []
        newseq = []
        for i in range(len(self.sequence)):
            if not self.sequence[i] == '-':
                newseq.append(self.sequence[i])
                idxs.append(i)
        return Sequence(newseq, self.alphabet, self.name, self.info, gappy = False), idxs

    def find(self, findme, gappy = False):
        """ Find the position of the specified symbol or sub-sequence """
        if gappy == False or self.gappy == False:
            return ''.join(self.sequence).find(findme)
        else: # if the sequence is gappy AND the function is called with gappy = True THEN run the find on the de-gapped sequence
            degapped, idxs = self.getDegapped()
            idx = ''.join(degapped).find(findme)
            return idxs[idx] if idx >= 0 else -1

    def connectFragment(self, newFragment):
        """

        :param newFragment: A sequence from the same contig as the original sequence
        :return: A combined sequence
        """
        if self.name.strip().strip('_qr') == newFragment.name.strip().strip('_qr'):
            # newFragment is before current fragment
            if newFragment.fragmentEnd < self.fragmentStart:
                ambiguousN = ''.join(['N']*(self.fragmentStart-newFragment.fragmentEnd-1))
                self.sequence = newFragment.sequence + ambiguousN + self.sequence
                self.fragmentStart = newFragment.fragmentStart
            # newFragment is after current fragment
            elif self.fragmentEnd < newFragment.fragmentStart:
                ambiguousN = ''.join(['N']*(newFragment.fragmentStart-self.fragmentEnd-1))
                self.sequence = self.sequence + ambiguousN + newFragment.sequence
                self.fragmentEnd = newFragment.fragmentEnd
            # newFragment is inside current sequence, so only replace ambiguous regions
            elif self.fragmentStart < newFragment.fragmentStart and self.fragmentEnd > newFragment.fragmentEnd:
                start = newFragment.fragmentStart - self.fragmentStart
                end = newFragment.fragmentEnd - self.fragmentStart
                subSeq = self.sequence[start:end+1]
                replaceSeq = ''
                replacement_cnt = 0

                # Replace any ambiguous characters in fragment region
                if len(subSeq) != len(newFragment.sequence):
                    print("differing sequence lengths?", len(subSeq), len(newFragment.sequence))
                for (i, sym) in enumerate(subSeq):
                    if sym == 'N':
                        replaceSeq += newFragment.sequence[i]
                        replacement_cnt += 1
                    else:
                        if self.percentIdentity >= newFragment.percentIdentity:
                            replaceSeq += sym
                        else:
                            replaceSeq += newFragment.sequence[i]
                            replacement_cnt += 1

                self.percentIdentity = (self.percentIdentity+newFragment.percentIdentity)/2
                # print("".join(self.sequence[:start]), replaceSeq, "".join(self.sequence[end+1:]))
                self.sequence = "".join(self.sequence[:start]) + replaceSeq + "".join(self.sequence[end+1:])
        else:
            print("Sequences are from different contigs:")
            print(self.name, newFragment.name)


class fastaFile():
    def __init__(self, fasta):
        """
        Create a vcfFile instance.
        :param entries: an iterable of entries or a filename
        """
        if isinstance(fasta, str): # filename
            entries = readFastaFile(fasta)
            self.contigs = dict()
            for entry in entries:
                self.contigs[entry.name] = entry
        else:
            self.contigs = dict()
            for entry in fasta:
                self.contigs[entry.name] = entry

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

"""
Below are some useful methods for loading data from strings and files.
Recognize the FASTA format (nothing fancy).
"""
def readFasta(string, alphabet = None, ignore = False, gappy = False, parse_defline = True):
    """ Read the given string as FASTA formatted data and return the list of
        sequences contained within it.
        If alphabet is specified, use it, if None (default) then guess it.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences.
        If parse_defline is False, the name will be set to everything before the first space, else parsing will be attempted."""
    seqlist = []    # list of sequences contained in the string
    seqname = None  # name of *current* sequence
    seqinfo = None
    seqdata = []    # sequence data for *current* sequence
    for line in string.splitlines():    # read every line
        if len(line) == 0:              # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence
            if seqname:     # check if we've got one current
                try:
                    current = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
                    seqlist.append(current)
                except RuntimeError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            # now collect data about the new sequence
            seqinfo = line[1:].split() # skip first char
            if len(seqinfo) > 0:
                try:
                    if parse_defline:
                        parsed = parseDefline(seqinfo[0])
                        seqname = parsed[0]
                        seqinfo = line[1:]
                    else: # we are not parsing the sequence name so no need to duplicate it in the info
                        seqname = seqinfo[0]
                        if len(seqinfo) > 0: # more than a name
                            edited_info = ''
                            for infopart in seqinfo[1:]:
                                edited_info += infopart + ' '
                            seqinfo = edited_info
                        else:
                            seqinfo = ''
                except IndexError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            else:
                seqname = ''
                seqinfo = ''
            seqdata = []
        else:               # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        try:
            lastseq = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
            seqlist.append(lastseq)
        except RuntimeError as errmsg:
            if not ignore:
                raise RuntimeError(errmsg)
    return seqlist

def parseDefline(string):
    """ Parse the FASTA defline (see http://en.wikipedia.org/wiki/FASTA_format)
        GenBank, EMBL, etc                gi|gi-number|gb|accession|locus
        SWISS-PROT, TrEMBL                sp|accession|name
        ...
        Return a tuple with
        [0] primary search key, e.g. UniProt accession, Genbank GI
        [1] secondary search key, e.g. UniProt name, Genbank accession
        [2] source, e.g. 'sp' (SwissProt/UniProt), 'tr' (TrEMBL), 'gb' (Genbank)
    """
    if len(string) == 0: return ('', '', '', '')
    s = string.split()[0]
    if re.match("^sp\|[A-Z][A-Z0-9]{5}\|\S+", s):            arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^tr\|[A-Z][A-Z0-9]*\|\S+", s): arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^gi\|[0-9]*\|\S+\|\S+", s):               arg = s.split('|');  return (arg[1], arg[3], arg[0], arg[2])
    elif re.match("gb\|\S+\|\S+", s):                        arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("emb\|\S+\|\S+", s):                       arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^refseq\|\S+\|\S+", s):                   arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    else: return (s, '', '', '')

def readFastaFile(filename, alphabet = None, ignore = False, gappy = False, parse_defline = True):
    """ Read the given FASTA formatted file and return the list of sequences
        contained within it. Note that if alphabet is NOT specified, it will take a
        separate guess for each sequence.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences.
        If parse_defline is False, the name will be set to everything before the first space, else parsing will be attempted."""

    seqlist = []
    with open(filename) as fh:
        seq = ''
        name = ''
        info = ''
        batch = '' # a batch of rows including one or more complete FASTA entries
        building = False
        for idx, row in enumerate(fh):
            row = row.strip()
            if len(row) > 0:
                if row.startswith('>'):
                    if building:
                        seqlist.append(Sequence(seq, name=name, info=info))
                    row = row.strip('>').split()
                    name = row[0]
                    info = " ".join(row)
                    seq = ''
                    building = True
                else:
                    seq += row
        seqlist.append(Sequence(seq, name=name, info=info))
    return seqlist

def writeFastaFile(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.writeFasta())
    fh.close()

def getMarkov(seqs, order = 0):
    """ Retrieve the Markov stats for a set of sequences. """
    myseqs = seqs
    if seqs is Sequence:
        myseqs = list([seqs])
    myalpha = None
    for seq in myseqs:
        if myalpha == None:
            myalpha = seq.alphabet
        else:
            if seq.alphabet != myalpha:
                raise RuntimeError('Sequence ' + seq.name + ' uses an invalid alphabet ')
    jp = Joint([myalpha for _ in range(order)])
    for seq in myseqs:
        for i in range(len(seq) - order):
            sub = seq[i:i + order + 1]
            jp.observe(sub)
    return jp

def getCount(seqs, findme = None):
    if findme != None:
        cnt = 0
        for seq in seqs:
            cnt += seq.count(findme)
        return cnt
    else:
        if len(seqs) > 0:
            alpha = seqs[0].alphabet
            patcnt = {}
            for a in alpha:
                patcnt[a] = getCount(seqs, a)
        return patcnt



if __name__ == '__main__':
    try:
        input_fasta = sys.argv[1]
        output_fasta = sys.argv[2]
        fasta = readFastaFile(input_fasta)
        seq_cnt = 1
        for seq in fasta:
            seq.info = ""
            seq.name = input_fasta.strip().split('/')[-1].split('.')[0]+"_"+str(seq_cnt)
            seq_cnt += 1

        writeFastaFile(output_fasta, fasta)
    except IndexError:
        print("Usage: <Input File> <Output File> \n",
              "Removes annotations from Fasta Headers \n")


