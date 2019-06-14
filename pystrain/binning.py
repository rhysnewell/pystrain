import pystrain.vcf as vcf
import pystrain.fastaindex as fastaindex
import pystrain.sequence as sequence
import sys, os, glob

def bin_contigs(vcf_file, fai_file):
    count_dict = {}
    for contig in vcf_file.contigs().keys():
        count_dict[contig] = vcf_file.getContigVariation(contig, fai_file.contigs[contig])