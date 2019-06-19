#!/usr/bin/env python3
import pystrain.coords as coords
import pystrain.sequence as sequence
import pyfaidx
import sys, os, glob
import importlib
import threading

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


# tests()

def binOnThread(assemblyCoords, contig, min_length, min_id, min_cov, min_genome_length, bin_cnt, simple=True, outputDirectory = "./"):
    spool = {}
    if assemblyCoords.query.index[contig].rlen >= min_length:

        for query_coords in assemblyCoords.generate(contig, source='q'):
            # Add in initial fragment if it is good
            if query_coords.percent_id >= min_id and query_coords.s2_len >= min_length and query_coords.q_cov >= min_cov:
                try:
                    if query_coords.seen:
                        continue
                except AttributeError:
                    query_coords.seen = False
                if simple:
                    spool[query_coords.q_tag] = assemblyCoords.query.fetch(query_coords.q_tag, 1,
                                                                           assemblyCoords.query.index[
                                                                               query_coords.q_tag].rlen)
                    query_coords.seen = True

                else:
                    try:
                        contig_fraction = assemblyCoords.query.fetch(query_coords.q_tag,
                                                                     query_coords.s2_start, query_coords.s2_end)
                        contig_fraction.percentIdentity = query_coords.percent_id
                        spool[query_coords.q_tag] = connectFragment(spool[query_coords.q_tag],
                                                                    contig_fraction)
                        query_coords.seen = True
                    except KeyError:
                        contig_fraction = assemblyCoords.query.fetch(query_coords.q_tag,
                                                                     query_coords.s2_start, query_coords.s2_end)
                        contig_fraction.percentIdentity = query_coords.percent_id
                        spool[query_coords.q_tag] = contig_fraction
                        query_coords.seen = True

            # Then begin search through reference contig
            source = 'r'
            tag = query_coords.r_tag
            searching = True
            matched = False
            while searching is True:
                for entry_coords in assemblyCoords.generate(tag, source=source):
                    if entry_coords.percent_id >= min_id and entry_coords.s2_len >= min_length and entry_coords.q_cov >= min_cov:
                        try:
                            if entry_coords.seen:
                                continue
                        except AttributeError:
                            entry_coords.seen = False

                        matched = True
                        # print(entry_coords)
                        if simple:
                            spool[entry_coords.q_tag] = assemblyCoords.query.fetch(entry_coords.q_tag, 1,
                                                                                   assemblyCoords.query.index[
                                                                                       entry_coords.q_tag].rlen)
                            entry_coords.seen = True

                        else:
                            try:
                                contig_fraction = assemblyCoords.query.fetch(entry_coords.q_tag,
                                                                             entry_coords.s2_start, entry_coords.s2_end)
                                contig_fraction.percentIdentity = entry_coords.percent_id
                                spool[entry_coords.q_tag] = connectFragment(spool[entry_coords.q_tag],
                                                                            contig_fraction)
                                entry_coords.seen = True
                            except KeyError:
                                contig_fraction = assemblyCoords.query.fetch(entry_coords.q_tag,
                                                                             entry_coords.s2_start, entry_coords.s2_end)
                                contig_fraction.percentIdentity = entry_coords.percent_id
                                spool[entry_coords.q_tag] = contig_fraction
                                entry_coords.seen = True
                        if source == 'r':
                            tag = entry_coords.q_tag
                            source = 'q'
                            break
                        else:
                            tag = entry_coords.r_tag
                            source = 'r'
                            break
                if matched is True:
                    matched = False
                else:
                    try:
                        if entry_coords.seen:
                            entry_coords.seen = False
                    except AttributeError:
                        entry_coords.seen = False
                    searching = False
        if len("".join([seq.seq for seq in spool.values()])) >= min_genome_length:
            bin_cnt += 1
            # bins["bin."+str(bin_cnt)] = spool
            print("Working on bin:", bin_cnt)
            with open(outputDirectory + "/" + str(bin_cnt) + ".fna", 'w') as fh:
                for contig in spool.keys():
                    seq = spool[contig]
                    fasta = ">" + seq.name + '\n'
                    fasta += seq.seq + '\n'
                    fh.write(fasta)

def binContigs(assemblyCoords, min_length=500, min_id=97, min_cov=5, min_genome_length=250000, simple=True, outputDirectory = "./", numThreads=5):
    bin_cnt = 0
    bins = {}
    threadLock = threading.Lock()
    threads = []
    for contig in assemblyCoords.query.index.keys():
        spool = {}
        if assemblyCoords.query.index[contig].rlen >= min_length:

            for query_coords in assemblyCoords.generate(contig, source='q'):
                # Add in initial fragment if it is good
                if query_coords.percent_id >= min_id and query_coords.s2_len >= min_length and query_coords.q_cov >= min_cov:
                    try:
                        if query_coords.seen:
                            continue
                    except AttributeError:
                        query_coords.seen = False
                    if simple:
                        spool[query_coords.q_tag] = assemblyCoords.query.fetch(query_coords.q_tag, 1,
                                                                                   assemblyCoords.query.index[query_coords.q_tag].rlen)
                        query_coords.seen = True

                    else:
                        try:
                            contig_fraction = assemblyCoords.query.fetch(query_coords.q_tag,
                                                                         query_coords.s2_start, query_coords.s2_end)
                            contig_fraction.percentIdentity = query_coords.percent_id
                            spool[query_coords.q_tag] = connectFragment(spool[query_coords.q_tag],
                                                                            contig_fraction)
                            query_coords.seen = True
                        except KeyError:
                            contig_fraction = assemblyCoords.query.fetch(query_coords.q_tag,
                                                                         query_coords.s2_start, query_coords.s2_end)
                            contig_fraction.percentIdentity = query_coords.percent_id
                            spool[query_coords.q_tag] = contig_fraction
                            query_coords.seen = True

                # Then begin search through reference contig
                source = 'r'
                tag = query_coords.r_tag
                searching = True
                matched = False
                while searching is True:
                    for entry_coords in assemblyCoords.generate(tag, source=source):
                        if entry_coords.percent_id >= min_id and entry_coords.s2_len >= min_length and entry_coords.q_cov >= min_cov:
                            try:
                                if entry_coords.seen:
                                    continue
                            except AttributeError:
                                entry_coords.seen = False

                            matched = True
                            # print(entry_coords)
                            if simple:
                                spool[entry_coords.q_tag] = assemblyCoords.query.fetch(entry_coords.q_tag, 1,
                                                                                       assemblyCoords.query.index[
                                                                                           entry_coords.q_tag].rlen)
                                entry_coords.seen = True

                            else:
                                try:
                                    contig_fraction = assemblyCoords.query.fetch(entry_coords.q_tag,
                                                                                 entry_coords.s2_start, entry_coords.s2_end)
                                    contig_fraction.percentIdentity = entry_coords.percent_id
                                    spool[entry_coords.q_tag] = connectFragment(spool[entry_coords.q_tag],
                                                                                    contig_fraction)
                                    entry_coords.seen = True
                                except KeyError:
                                    contig_fraction = assemblyCoords.query.fetch(entry_coords.q_tag,
                                                                                 entry_coords.s2_start, entry_coords.s2_end)
                                    contig_fraction.percentIdentity = entry_coords.percent_id
                                    spool[entry_coords.q_tag] = contig_fraction
                                    entry_coords.seen = True
                            if source == 'r':
                                tag = entry_coords.q_tag
                                source = 'q'
                                break
                            else:
                                tag = entry_coords.r_tag
                                source = 'r'
                                break
                    if matched is True:
                        matched = False
                    else:
                        try:
                            if entry_coords.seen:
                                entry_coords.seen = False
                        except AttributeError:
                            entry_coords.seen = False
                        searching = False
            if len("".join([seq.seq for seq in spool.values()])) >= min_genome_length:
                bin_cnt += 1
                # bins["bin."+str(bin_cnt)] = spool
                print("Working on bin:", bin_cnt)
                with open(outputDirectory + "/" + "bin." + str(bin_cnt) + ".fna", 'w') as fh:
                    for contig in spool.keys():
                        seq = spool[contig]
                        fasta = ">" + seq.name + '\n'
                        fasta += seq.seq + '\n'
                        fh.write(fasta)
    return bins


def buildContigs(assemblyCoords, oldBin, simple=True, outputDirectory='./', min_length = 2000, min_id = 85, min_cov=90):
    """

    :param assemblyCoords: The alignment coords produced by nucmer when aligning two assemblies contigs together
    :param binCoords:  a metagenomic bin
    :return: Two new bin .fna files with any additional contigs that were found to be present in both samples and
    binned in one species
    """
    new_contigs = {}
    connections = {}
    for contig in oldBin.index.keys():
        for query_contig in assemblyCoords.generate(contig, source='q'):
            for ref_contig in assemblyCoords.generate(query_contig.r_tag, source='r'):
                if ref_contig.q_tag not in oldBin.index.keys():
                    dist = ref_contig.s1_start - query_contig.s1_end
                    if simple:
                        if ref_contig.percent_id >= min_id and ref_contig.s2_len >= min_length and ref_contig.q_cov >= min_cov:
                            new_contigs[ref_contig.q_tag] = assemblyCoords.query.fetch(ref_contig.q_tag, 1,
                                                                                       assemblyCoords.query.index[ref_contig.q_tag].rlen)
                    else:
                        if ref_contig.percent_id >= min_id and ref_contig.s2_len >= min_length and ref_contig.q_cov >= min_cov:
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

    return new_contigs


def twoSampleBuildContigs(assemblyCoordsFile, binCoordsDirectory, outputDirectory, minLength, minMatch, minCov, simple):
    binCoords = glob.glob(binCoordsDirectory+"/*")
    assembly_coords = coords.readCoordFile(assemblyCoordsFile)

    try:
        os.mkdir(outputDirectory)
    except FileExistsError:
        print("Overwriting existing files")

    for file in binCoords:
        if file.endswith(".fna") or file.endswith(".fa"):
            print("Working on: ", file)
            oldBin = pyfaidx.Faidx(file)
            new_contigs = buildContigs(assembly_coords, oldBin, simple, outputDirectory, minLength, minMatch, minCov)
            with open(outputDirectory + "/new_" + file.split('/')[-1], 'w') as fh:

                for contig in oldBin.index.keys():
                    seq = oldBin.fetch(contig, 1, oldBin.index[contig].rlen)
                    fasta = ">" + seq.name + '\n'
                    fasta += seq.seq + '\n'
                    fh.write(fasta)
                for contig in new_contigs.keys():
                    fasta = '>' + contig + '\n'
                    fasta += new_contigs[contig].seq
                    fasta += '\n'
                    fh.write(fasta)

    print("done!")

def spooledBinning(assemblyCoordFile, outputDirectory,  minLength, minMatch, minCov):
    assembly_coords = coords.readCoordFile(assemblyCoordFile)

    try:
        os.mkdir(outputDirectory)
    except FileExistsError:
        print("Overwriting existing files")

    binContigs(assembly_coords, minLength, minMatch, minCov, outputDirectory=outputDirectory)
    # for mag in bins.keys():
    #     print("Working on: ", mag)
    #     with open(outputDirectory + "/" + mag+".fna", 'w') as fh:
    #         for contig in bins[mag].keys():
    #             seq = bins[mag][contig]
    #             fasta = ">" + seq.name + '\n'
    #             fasta += seq.seq + '\n'
    #             fh.write(fasta)
    print("done!")


if __name__ == "__main__":
    try:
        mode = sys.argv[1]
        if mode == 'build':
            try:
                assembly = sys.argv[2]
                bins = sys.argv[3]
                output_directory = sys.argv[4]
                min_length = sys.argv[5]
                min_match = sys.argv[6]
                min_query_cov = sys.argv[7]
                try:
                    present = sys.argv[8]
                    simple = False
                except IndexError:
                    simple = True
                twoSampleBuildContigs(assembly, bins, output_directory, float(min_length), float(min_match),
                                      float(min_query_cov), simple)

            except IndexError:
                print(
                    "Usage: completecontigs.py build <AssemblyCoords> <BinsDirectory> <OutputDirectory> <MinimumMatchLength> <MinimumMatchID> <MinimumQueryCoverage> [ComplexMode]")

        if mode == 'bin':
            try:
                assembly = sys.argv[2]
                output_directory = sys.argv[3]
                min_length = sys.argv[4]
                min_match = sys.argv[5]
                min_query_cov = sys.argv[6]
                spooledBinning(assembly, output_directory, float(min_length), float(min_match),
                                      float(min_query_cov))

            except IndexError:
                print(
                    "Usage: completecontigs.py bin <AssemblyCoords> <OutputDirectory> <MinimumMatchLength> <MinimumMatchID> <MinimumQueryCoverage>")
        else:
            print("Usage: completecontigs.py <Mode> \n"
                  "Mode: bin OR build")


    except IndexError:
        print("Usage: completecontigs.py <Mode> \n"
              "Mode: bin OR build")
