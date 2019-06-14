import shlex
import pystrain.ival as ival
import pystrain.fastaindex as fastaindex

class vcfEntry():

    '''
    GFF fields:
    seqname - The name of the sequence. Must be a chromosome or scaffold.
    source - The program that generated this feature.
    feature - The name of this type of feature. Some examples of standard feature types are "CDS" "start_codon" "stop_codon" and "exon"li>
    start - The starting position of the feature in the sequence. The first base is numbered 1.
    end - The ending position of the feature (inclusive).
    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ":.":.
    strand - Valid entries include "+", "-", or "." (for don't know/don't care).
    frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be ".".
    group - All lines with the same group are linked together into a single item.
    '''
    def __init__(self, contig, pos, ID, ref, alt, qual, filters, info_dict):
        self.contig = contig
        self.pos = pos
        self.ID = ID
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filters = filters
        self.info_dict = info_dict
        # self.formats = formats
        # self.unknown = unknown

    def __getitem__(self, item):
        return self.info_dict[item]

    def __contains__(self, item):
        return item in self.info_dict




class vcfFile:
    """ Read GTF/GFF file.

        See http://genome.ucsc.edu/FAQ/FAQformat#format1
    """

    def __init__(self, entries, filter_feature = None):
        """
        Create a vcfFile instance.
        :param entries: an iterable of entries or a filename
        """
        if isinstance(entries, str): # filename
            self.contigs = readvcfFile(entries, filter_feature)
        else:
            self.contigs = dict()
            for entry in entries:
                # check if the chromosome has been seen before
                try:
                    self.contigs[entry.contig].append(entry)
                except KeyError:
                    self.contigs[entry.contig] = [entry]
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
        if isinstance(item, vcfEntry):
            tree = self.contigs.get(item.contig)
            if tree == None: return False
            else: return True
        else:
            return False

    def getContigVariation(self, contig, contig_idx):
        contig_entries = self.generate(contig)
        variation_count = 0
        for entry in contig_entries:
            variation_count += 1

        return variation_count/contig_idx.length

    def getIndex(self, group_attribute):
        """
        Create a dictionary for a specified attribute in the group list, e.g. "gene_id", "gene_name" or "transcript_id"
        :param group_attribute:
        :return: a dictionary keyed by the values of the nominated group_attribute,
        pointing to the entry in the chromosome interval-tree
        """
        d = {}
        for entry in self:
            if group_attribute in entry.info_dict: # is the entry indexed with the attribute
                if not entry.info_dict[group_attribute] in d: # only add the first entry for this particular value
                    d[entry.info_dict[group_attribute]] = entry
        return d

def readvcfFile(filename, filter_feature = None):
    """ Read a GTF/GFF/VCF file.
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
                continue # ignore empty lines
            if words[0].strip().startswith('#'):
                continue # comment
            if words[0].strip().startswith('browser'):
                continue # ignore
            if words[0].strip().startswith('track'):
                continue # ignore
            try:
                seqname = words[0]
                pos = words[1]
                ID = words[2]
    
                ref = words[3]
                alt = words[4]
                qual = words[5]
    
                filters = words[6]
                info = words[7]
                info_dict = {}
                for item in info.split(';'):
                    i = item.split('=')

                    try:
                        info_dict[i[0]] = i[1]
                    except IndexError:
                        info_dict[i[0]] = i[0]

                for type, value in zip(words[8].split(':'), words[9].split(':')):
                    info_dict[type] = value

                # formats = words[8]
                # unknown = words[9]
    
                entry = vcfEntry(seqname, pos, ID, ref, alt, qual, filters, info_dict)
                # check if the chromosome has been seen before
                try:
                    contigs[seqname].append(entry)
                except KeyError:
                    contigs[seqname] = [entry]
                
            except RuntimeError as e:
                if not acceptHeaderRows:
                    raise RuntimeError('Error in VCF file at row %d (%s)' % (row, e))
                else:
                    headerRow = words
                    acceptHeaderRows -= 1 # count down the number of header rows that can occur

    return contigs

def writevcfFile(entries, filename, header = None):
    """ Save the GTF entries to a file.
    """
    f = open(filename, 'w')
    if header:
        f.write(header + '\n')
    for row in entries:
        f.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s" % (row.contig, row.source, row.feature, row.start, row.end, row.score, row.strand, row.frame, row.group))
        f.write("\n")
    f.close()

if __name__ == '__main__':
    bf = vcfFile('tests/r2.parent.d182.vcf')
    fai = fastaindex.faiFile('tests/r2.parent.d182.assembly.fna.fai')
    print(bf.contigs.keys())
    g = bf.generate('k141_1')
    print(next(g))
    print(next(g))
    print(next(g))
    cnt = 0
    print(bf.getContigVariation('k141_1', fai.contigs['k141_1']))
    for contig in bf.contigs.keys():
        print(bf.getContigVariation(contig, fai.contigs[contig]))
