import pystrain.vcf as vcf
import pystrain.fastaindex as fastaindex
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import pystrain.sequence as sequence
import sys, os, glob
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfaidx
# import pysam
import collections
import scipy.cluster.hierarchy as hcluster


#Keeps exiting with error code 139 whenever iterated through?
#Possible bug in pysam
# bam = pysam.AlignmentFile("tests/7seqs.reads_for_seq1_and_seq2.bam", "rb")
# index = pysam.IndexedReads(bam)
# contigs = bam.references
#
# pileups = bam.pileup()
# variations = {}
# for p in pileups:
#     print("\n"+p.reference_name+"\t"+p.reference_name+"\t"+ str(p.reference_pos))
#
#     queries = p.get_query_sequences(mark_matches=True, add_indels=True)
    # print("\ncoverage at base %s = %s" %
    #       (p.pos, p.n))
    # for pileupread in p.pileups:
    #     if not pileupread.is_del and not pileupread.is_refskip:
    #         query position is None if is_del or is_refskip is set.
        # print('\tbase in read %s = %s' %
        #       (pileupread.alignment.query_name,
        #        pileupread.alignment.query_sequence[pileupread.query_position]))
        # print(pileupread.alignment.get_reference_sequence()[pileupread.query_position])

# bf = vcf.vcfFile('tests/paludibacter.vcf')
# fai = fastaindex.faiFile('tests/paludibacter.binned_contigs.38.fa.fai')

bf = vcf.vcfFile('tests/r2.parent.d182.vcf')
fai = fastaindex.faiFile('tests/r2.parent.d182.assembly.fna.fai')

# variants = pyfaidx.FastaVariant('tests/r2.parent.d182.assembly.fna', 'tests/r2.parent.d182.vcf', het=False, hom=False)

def bin_contigs(vcf_file, fai_file):
    value_dict = collections.OrderedDict()
    vc = []
    lengths = []
    for contig in vcf_file.contigs.keys():
        vpfl = vcf_file.getContigVariation(contig, fai_file.contigs[contig], 100)  # Variation Per Fragment Length
        value_dict[contig] = [np.std(vpfl), np.mean(vpfl), np.median(vpfl)]

    # return pd.DataFrame.from_dict(value_dict, orient='index')
    return value_dict

df_v = bin_contigs(bf, fai)

values = np.empty([len(df_v.keys()), 3])

count = 0
for (item, value) in df_v.items():
    values[count] = value
    count += 1

values[:,0] = (values[:,0]-min(values[:,0]))/(max(values[:,0] - min(values[:,0])))
values[:,1] = (values[:,1]-min(values[:,1]))/(max(values[:,1] - min(values[:,1])))
values[:,2] = (values[:,2]-min(values[:,2]))/(max(values[:,2] - min(values[:,2])))

# clustering = AgglomerativeClustering(n_clusters=None, compute_full_tree=True, distance_threshold=0.01)
# clustering.fit_predict(values[0:100])


clusters = hcluster.fclusterdata(values[0:40000], 0.01, criterion="distance")
plt.scatter(*np.transpose(values[0:40000]), c=clusters)
plt.axis("equal")
title = "threshold: %f, number of clusters: %d" % (0.01, len(set(clusters)))
plt.title(title)
plt.show()
#
# kmeans_model = KMeans(n_clusters=20, random_state=1).fit(df_v.iloc[:, :])
# labels = kmeans_model.predict(df_v)
# scatter = plt.scatter(df_v[1], df_v[2], c=kmeans_model.labels_)
# plt.show()
#
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
#
# Z = linkage(df_v, 'ward')
# fig = plt.figure()
# dn = dendrogram(Z)
# plt.show()
# fig.savefig("linkagebin_cluster.pdf", bbox_inches="tight")
