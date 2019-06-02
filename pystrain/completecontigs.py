import pystrain.coords as coords
import pystrain.sequence as sequence
import importlib
importlib.reload(sequence)
importlib.reload(coords)


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

assembly_coords = coords.readCoordFile("tests/test_assemblies.coords")
bin_coords = coords.readCoordFile("tests/test_bins.coords")

def buildContigs(assemblyCoords, binCoords):
    # tag_dict = {}
    contig_occ = {}
    for contig in assemblyCoords.contigs:
        if assemblyCoords.source[contig] == 'r':
            for entry in assemblyCoords.generate(contig):
                if entry.q_tag in binCoords.contigs.keys() and entry.r_tag not in binCoords.contigs.keys():
                    contig_fraction = assemblyCoords.reference[entry.r_tag][entry.s1_start:entry.s1_end]
                    binCoords.reference[entry.r_tag] = sequence.Sequence(contig_fraction, name=entry.r_tag)
                    print(contig_fraction)
    sequence.writeFastaFile('ref_new_bin.fna', binCoords.reference.values())


buildContigs(assembly_coords, bin_coords)


