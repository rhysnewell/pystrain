#!/usr/bin/env python3
import pystrain.coords as coords
import pystrain.sequence as sequence
import sys, os, glob
import importlib

def tests():
    seq_a = sequence.Sequence('ATGCGT', name = 'contig1', fragmentposition=[0, 5])
    print(seq_a, seq_a.fragmentStart, seq_a.fragmentEnd)
    seq_b = sequence.Sequence('GCGAGT', name = 'contig1', fragmentposition=[12, 17])
    seq_a.connectFragment(seq_b)
    # Connect fragment after original
    print(seq_a, seq_a.fragmentStart, seq_a.fragmentEnd)
    # Wrong contig
    seq_c = sequence.Sequence('CCCCCC', name = 'cotnig1', fragmentposition=[6, 11])
    seq_a.connectFragment(seq_c)
    seq_c = sequence.Sequence('CCCCCC', name = 'contig1', fragmentposition=[6, 11])
    seq_a.connectFragment(seq_c)
    # Substitute fragment in centre
    print(seq_a, seq_a.fragmentStart, seq_a.fragmentEnd)
    seq_a = sequence.Sequence('ATGCGT', name='contig1', fragmentposition=[0, 5])
    seq_b = sequence.Sequence('GCGAGT', name='contig1', fragmentposition=[12, 17])
    seq_b.connectFragment(seq_a)
    # Connect fragment before original
    print(seq_b, seq_b.fragmentStart, seq_b.fragmentEnd)

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


def buildContigs(assemblyCoords, binCoords, simple=True, outputDirectory='./'):
    """

    :param assemblyCoords: The alignment coords produced by nucmer when aligning two assemblies contigs together
    :param binCoords:  As above, but two bins contigs
    :return: Two new bin .fna files with any additional contigs that were found to be present in both samples and
    binned in one species
    """
    new_ref_contigs = {}
    new_query_contigs = {}
    for contig in assemblyCoords.contigs:
        if assemblyCoords.source[contig] == 'r':
            for entry in assemblyCoords.generate(contig):
                if entry.q_tag in binCoords.query.keys() and entry.r_tag not in binCoords.reference.keys():
                    contig_fraction = assemblyCoords.reference[entry.r_tag][entry.s1_start:entry.s1_end+1]
                    if simple:
                        new_ref_contigs[entry.r_tag] = assemblyCoords.reference[entry.r_tag]
                    else:
                        try:
                            new_ref_contigs[entry.r_tag].connectFragment(
                                sequence.Sequence(contig_fraction,
                                                  name=entry.r_tag,
                                                  fragmentposition = [entry.s1_start, entry.s1_end]))
                        except KeyError:
                            new_ref_contigs[entry.r_tag] = sequence.Sequence(contig_fraction,
                                                                             name=entry.r_tag,
                                                                             fragmentposition = [entry.s1_start, entry.s1_end])
        if assemblyCoords.source[contig] == 'q':
            for entry in assemblyCoords.generate(contig):
                if entry.q_tag not in binCoords.query.keys() and entry.r_tag in binCoords.reference.keys():
                    contig_fraction = assemblyCoords.query[entry.q_tag][entry.s2_start:entry.s2_end+1]
                    if simple:
                        new_query_contigs[entry.q_tag] = assemblyCoords.query[entry.q_tag]
                    else:
                        try:
                            new_query_contigs[entry.q_tag].connectFragment(
                                sequence.Sequence(contig_fraction,
                                                  name=entry.q_tag,
                                                  fragmentposition = [entry.s2_start, entry.s2_end]))
                        except KeyError:
                            new_query_contigs[entry.q_tag] = sequence.Sequence(contig_fraction,
                                                                             name=entry.q_tag,
                                                                             fragmentposition = [entry.s2_start, entry.s2_end])
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
        print("Overwriting previous query bins")

    sequence.writeFastaFile(outputDirectory+"/referencebins/"+"new_"+binCoords.reference_name,
                            list(binCoords.reference.values())+list(new_ref_contigs.values()))

    sequence.writeFastaFile(outputDirectory+"/querybins/"+"new_"+binCoords.query_name,
                            list(binCoords.query.values())+list(new_query_contigs.values()))



# buildContigs(assembly_coords, bin_coords)

def twoSampleBuildContigs(assemblyCoordsFile, binCoordsDirectory, outputDirectory):
    binCoords = glob.glob(binCoordsDirectory+"/*.coords")
    assembly_coords = coords.readCoordFile(assemblyCoordsFile)
    for file in binCoords:
        if file.endswith(".coords"):
            bin_coords = coords.readCoordFile(file)
            buildContigs(assembly_coords, bin_coords, True, outputDirectory)

if __name__ == "__main__":
    try:
        assembly = sys.argv[1]
        bins = sys.argv[2]
        try:
            output_directory = sys.argv[3]
        except IndexError:
            output_directory = './'
        twoSampleBuildContigs(assembly, bins, output_directory)
    except IndexError:
        print("Usage: completecontigs.py <AssemblyCoords> <BinCoordsDirectory> <OutputDirectory>")
