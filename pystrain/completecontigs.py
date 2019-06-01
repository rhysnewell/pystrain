import pystrain.coords as coords
import pystrain.sequence as sequence
import importlib
importlib.reload(sequence)
importlib.reload(coords)


testcoords = coords.readCoordFile("tests/out.coords")
ref = sequence.readFastaFile("tests/73.20100900_E3D.15.fna")
query = sequence.readFastaFile("tests/r2.parent.d77.174_unicyc.fna")

ref_dict = {}
for seq in ref:
    ref_dict[seq.name] = seq

query_dict = {}
for seq in query:
    query_dict[seq.name] = seq

def buildContigs(contig_coords):
    # tag_dict = {}
    new_seq_dict = {}
    for contig in contig_coords.contigs:
        # if contig_coords.source[contig] == 'r':
        #     print(contig)
        found_in = 0
        contains = 0
        for entry in contig_coords.generate(contig):
            print(entry)
            if entry.r_len > entry.q_len:
                contains += 1
            else:
                found_in += 1
        break

buildContigs(testcoords)


