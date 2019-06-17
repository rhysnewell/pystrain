import pystrain.vcf as vcf
import pystrain.fastaindex as fastaindex
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import pystrain.sequence as sequence
import sys, os, glob
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

bf = vcf.vcfFile('tests/r2.parent.d182.vcf')
fai = fastaindex.faiFile('tests/r2.parent.d182.assembly.fna.fai')

def bin_contigs(vcf_file, fai_file):
    value_dict = {}
    vc = []
    lengths = []
    for contig in vcf_file.contigs.keys():
        vpfl = vcf_file.getContigVariation(contig, fai_file.contigs[contig])  # Variation Per Fragment Length
        value_dict[contig] = [np.std(vpfl), np.mean(vpfl), np.median(vpfl)]

    return pd.DataFrame.from_dict(value_dict, orient='index')

df_v = bin_contigs(bf, fai)

# ks = []
# inertias = []
# for k in range(1, 20):
#     # Create a kmeans model on our data, using k clusters.  random_state helps ensure that the algorithm returns the same results each time.
#     kmeans_model = KMeans(n_clusters=k, random_state=1).fit(df_v.iloc[:, :])
#
#     # These are our fitted labels for clusters -- the first cluster has label 0, and the second has label 1.
#     labels = kmeans_model.labels_
#
#     # Sum of distances of samples to their closest cluster center
#     interia = kmeans_model.inertia_
#     inertias.append(interia)
#     ks.append(k)
#     print("k:", k, " cost:", interia)
#
# plt.plot(ks, inertias, 'bo')
# plt.show()

Z = linkage(df_v, 'ward')
fig = plt.figure()
dn = dendrogram(Z)
plt.show()
fig.savefig("linkagebin_cluster.pdf", bbox_inches="tight")
