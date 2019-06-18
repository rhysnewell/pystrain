#!/usr/bin/env python3
import pystrain.coords as coords
import pystrain.sequence as sequence
import pyfaidx
import sys, os, glob
import importlib

def connectFragment(oldFragment, newFragment):
    """

    :param newFragment: A sequence from the same contig as the original sequence from sequence class in pyfaidx
    :return: A combined sequence
    """

    # Set percentIdentity attribute placehodler in case it doesn't exist yet
    try:
        oldFragment.percentIdentity
    except AttributeError:
        oldFragment.percentIdentity = 100

    try:
        newFragment.percentIdentity
    except AttributeError:
        newFragment.percentIdentity = 100

    if oldFragment.name.strip() == newFragment.name.strip():
        # newFragment is before current fragment
        if newFragment.end < oldFragment.start:
            ambiguousN = ''.join(['N']*(oldFragment.start-newFragment.end-1))
            oldFragment.seq = newFragment.seq + ambiguousN + oldFragment.seq
            oldFragment.start = newFragment.start
        # newFragment is after current fragment
        elif oldFragment.end < newFragment.start:
            ambiguousN = ''.join(['N']*(newFragment.start-oldFragment.end-1))
            oldFragment.seq = oldFragment.seq + ambiguousN + newFragment.seq
            oldFragment.end = newFragment.end
        # newFragment is inside current sequence, so only replace ambiguous regions
        elif oldFragment.start < newFragment.start and oldFragment.end > newFragment.end:
            start = newFragment.start - oldFragment.start
            end = newFragment.end - oldFragment.start
            subSeq = oldFragment.seq[start:end+1]
            replaceSeq = ''
            replacement_cnt = 0

            # Replace any ambiguous characters in fragment region
            if len(subSeq) != len(newFragment.seq):
                print("differing sequence lengths?", len(subSeq), len(newFragment.seq))
            for (i, sym) in enumerate(subSeq):
                if sym == 'N':
                    replaceSeq += newFragment.seq[i]
                    replacement_cnt += 1
                else:

                    if oldFragment.percentIdentity >= newFragment.percentIdentity:
                        replaceSeq += sym
                    else:
                        replaceSeq += newFragment.seq[i]
                        replacement_cnt += 1

            oldFragment.percentIdentity = (oldFragment.percentIdentity+newFragment.percentIdentity)/2
            # print("".join(oldFragment.seq[:start]), replaceSeq, "".join(oldFragment.seq[end+1:]))
            oldFragment.seq = "".join(oldFragment.seq[:start]) + replaceSeq + "".join(oldFragment.seq[end+1:])
        return oldFragment
    else:
        print("Sequences are from different contigs:")
        print(oldFragment.name, newFragment.name)
        return oldFragment

def tests():
    seq_a = pyfaidx.Sequence('contig1', 'ATGCGT', start=1, end=6)
    print(seq_a, seq_a.start, seq_a.end)
    seq_b = pyfaidx.Sequence('contig1', 'GCGAGT', start=13, end=18)
    seq_a = connectFragment(seq_a, seq_b)
    # Connect fragment after original
    print(seq_a, seq_a.start, seq_a.end)
    # Wrong contig
    seq_c = pyfaidx.Sequence('different_contig1', 'CCCCCC', start=7, end=12)
    seq_a = connectFragment(seq_a, seq_c)
    seq_c = pyfaidx.Sequence('contig1', 'CCCCCC', start=7, end=12)
    seq_a = connectFragment(seq_a, seq_c)
    # Substitute fragment in centre
    print(seq_a, seq_a.start, seq_a.end)
    seq_a = pyfaidx.Sequence('contig1', 'ATGCGT', start=1, end=6)
    seq_b = pyfaidx.Sequence('contig1', 'GCGAGT', start=13, end=18)
    seq_b.percentIdentity = 97.1
    seq_b = connectFragment(seq_b, seq_a)
    # Connect fragment before original
    print(seq_b, seq_b.start, seq_b.end, seq_b.percentIdentity)
    seq_c = pyfaidx.Sequence('contig1', 'CCC', start=5, end=7)
    seq_c.percentIdentity = 97.2
    seq_b = connectFragment(seq_b, seq_c)
    print(seq_b, seq_b.start, seq_b.end, seq_b.percentIdentity)
    seq_c = pyfaidx.Sequence('contig1', 'AATTT', start=8, end=12)
    seq_b = connectFragment(seq_b, seq_c)
    print(seq_b, seq_b.start, seq_b.end, len(seq_b))


tests()
# testcoords = coords.readCoordFile("tests/out.coords")
# ref = sequence.readFastaFile("tests/73.20100900_E3D.15.fna")
# query = sequence.readFastaFile("tests/r2.parent.d77.174_unicyc.fna")

#
#
# ref_dict = {}
# for seq in ref:
#     ref_dict[seq.name] = seq
#
# query_dict = {}
# for seq in query:
#     query_dict[seq.name] = seq

# assembly_coords = coords.readCoordFile("tests/smallassembly.filter.coords")
# bin_coords = coords.readCoordFile("tests/bin.filter.coords")


def buildContigs(assemblyCoords, bin, simple=True, outputDirectory='./', min_length = 2000, min_id = 85):
    """

    :param assemblyCoords: The alignment coords produced by nucmer when aligning two assemblies contigs together
    :param binCoords:  a metagenomic bin
    :return: Two new bin .fna files with any additional contigs that were found to be present in both samples and
    binned in one species
    """
    new_contigs = {}
    connections = {}
    for contig in bin.index.keys():
        for query_contig in assemblyCoords.generate(contig, source='q'):
            for ref_contig in assemblyCoords.generate(query_contig.r_tag, source='r'):
                if ref_contig.q_tag not in bin.index.keys():
                    if simple:
                        if ref_contig.percent_id >= min_id and ref_contig.s2_len >= min_length:
                            new_contigs[ref_contig.q_tag] = assemblyCoords.query.fetch(ref_contig.q_tag, 1,
                                                                                       assemblyCoords.query.index[ref_contig.q_tag].rlen)
                    else:
                        if ref_contig.percent_id >= min_id and ref_contig.s2_len >= min_length:
                            try:
                                contig_fraction = assemblyCoords.query.fetch(ref_contig.q_tag,
                                                                             ref_contig.s2_start, ref_contig.s2_end)
                                contig_fraction.percentIdentity = ref_contig.percent_id
                                new_contigs[ref_contig.q_tag] = connectFragment(new_contigs[ref_contig.q_tag],
                                                                                contig_fraction)
                            except KeyError:
                                contig_fraction = assemblyCoords.query.fetch(ref_contig.q_tag,
                                                                             ref_contig.s2_start, ref_contig.s2_end)
                                contig_fraction.percentIdentity = ref_contig.percent_id
                                new_contigs[ref_contig.q_tag] = contig_fraction
        # if assemblyCoords.source[contig] == 'q':
        #     for entry in assemblyCoords.generate(contig):
        #         if entry.q_tag not in bin.query.keys() and entry.r_tag in bin.reference.keys():
        #
        #             if simple:
        #                 if entry.percent_id >= min_id and entry.s2_len >= min_length:
        #                     new_query_contigs[entry.q_tag] = assemblyCoords.query[entry.q_tag]
        #             else:
        #                 if entry.percent_id >= min_id and entry.s2_len >= min_length:
        #                     try:
        #                         new_query_contigs[entry.q_tag].connectFragment(
        #                             sequence.Sequence(contig_fraction,
        #                                               name=entry.q_tag,
        #                                               fragmentposition=[entry.s2_start, entry.s2_end],
        #                                               percentid=entry.percent_id))
        #                     except KeyError:
        #                         contig_fraction = assemblyCoords.query[entry.q_tag][entry.s2_start:entry.s2_end + 1]
        #                         new_query_contigs[entry.q_tag] = sequence.Sequence(contig_fraction,
        #                                                                            name=entry.q_tag,
        #                                                                            fragmentposition = [entry.s2_start, entry.s2_end],
        #                                                                            percentid=entry.percent_id)


    return new_contigs


def twoSampleBuildContigs(assemblyCoordsFile, binCoordsDirectory, outputDirectory, minLength, minMatch, simple):
    binCoords = glob.glob(binCoordsDirectory+"/*")
    assembly_coords = coords.readCoordFile(assemblyCoordsFile)

    try:
        os.mkdir(outputDirectory)
    except FileExistsError:
        print("Overwriting existing files")
    try:
        os.mkdir(outputDirectory+"/querybins")
    except FileExistsError:
        print("Overwriting previous query bins")
    try:
        os.mkdir(outputDirectory + "/referencebins")
    except FileExistsError:
        print("Overwriting previous reference bins")

    for file in binCoords:
        if file.endswith(".fna") or file.endswith(".fa"):
            print("Working on: ", file)
            bin = pyfaidx.Faidx(file)
            new_contigs = buildContigs(assembly_coords, bin, simple, outputDirectory, minLength, minMatch)
            with open(outputDirectory + "/new_" + file.split('/')[-1], 'w') as fh:

                for contig in bin.index.keys():
                    seq = bin.fetch(contig, 1, bin.index[contig].rlen)
                    fasta = ">" + seq.name + '\n'
                    fasta += seq.seq + '\n'
                    fh.write(fasta)
                for contig in new_contigs.keys():
                    fasta = '>' + contig + '\n'
                    fasta += new_contigs[contig].seq
                    fasta += '\n'
                    fh.write(fasta)

    print("done!")


if __name__ == "__main__":
    try:
        assembly = sys.argv[1]
        bins = sys.argv[2]
        output_directory = sys.argv[3]
        min_length = sys.argv[4]
        min_match = sys.argv[5]
        try:
            present = sys.argv[6]
            simple = False
        except IndexError:
            simple = True

        twoSampleBuildContigs(assembly, bins, output_directory, float(min_length), float(min_match), simple)
    except IndexError:
        print("Usage: completecontigs.py <AssemblyCoords> <BinsDirectory> <OutputDirectory> <MinimumMatchLength> <MinimumMatchID> [ComplexMode]")
