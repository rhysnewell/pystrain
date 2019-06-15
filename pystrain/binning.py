import pystrain.vcf as vcf
import pystrain.fastaindex as fastaindex
import pystrain.sequence as sequence
import sys, os, glob
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

bf = vcf.vcfFile('tests/r2.parent.d182.vcf')
fai = fastaindex.faiFile('tests/r2.parent.d182.assembly.fna.fai')

def bin_contigs(vcf_file, fai_file):
    count_dict = {}
    vc = []
    lengths = []
    for contig in vcf_file.contigs.keys():
        var = vcf_file.getContigVariation(contig, fai_file.contigs[contig])
        print(var, fai_file.contigs[contig].length, len(var), contig)
        # length = fai_file.contigs[contig].length
        # vc.append(var)
        # lengths.append(length)
        count_dict[contig] = var

    # sorted_values = sorted(count_dict.values())
    # # print(sorted_values)
    # plt.plot(sorted_values, 'bo')
    # plt.show()

bin_contigs(bf, fai)