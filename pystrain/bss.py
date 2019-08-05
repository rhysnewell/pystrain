import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib import cm
from scipy import signal
from pylab import *
from skimage import data, io, color
# from skimage.viewer import ImageViewer
from collections import OrderedDict
import nimfa
import pyfaidx

from sklearn.decomposition import FastICA, PCA, NMF
# #############################################################################
# Plot results

def plot_samples(S, axis_list=None):
    plt.scatter(S[:, 0], S[:, 1], s=2, marker='o', zorder=10,
                color='steelblue', alpha=0.5)
    if axis_list is not None:
        colors = ['orange', 'red']
        for color, axis in zip(colors, axis_list):
            axis /= axis.std()
            x_axis, y_axis = axis
            # Trick to get legend to work
            plt.plot(0.1 * x_axis, 0.1 * y_axis, linewidth=2, color=color)
            plt.quiver(0, 0, x_axis, y_axis, zorder=11, width=0.01, scale=6,
                       color=color)

    plt.hlines(0, -3, 3)
    plt.vlines(0, -3, 3)
    plt.xlim(-3, 3)
    plt.ylim(-3, 3)
    plt.xlabel('x')
    plt.ylabel('y')


def melt_dict(strain_dict):
    array = np.zeros(shape=(0,0))
    for key in sorted(strain_dict.keys()):
        # print(key, strain_dict[key])
        for bases in zip(*strain_dict[key]):
            print(bases)
            array = np.append(array, bases)
    return array


def read_strainm(filenames):
    # array = [[] for i in filenames]
    array = {}
    for (idx, filename) in enumerate(filenames):
        with open(filename) as f:
            for line in f:
                if line.startswith('pos'):
                    continue
                else:
                    line = line.strip().split()
                    # hexdec = int('a'.join(line[1:5]) + 'a', 16)
                    # depth = int(line[7])
                    # abundances = [int(x)/depth for x in line[1:5]]
                    try:
                        # multi file
                        # array[int(line[0])][i] = int(line[7])
                        # single file
                        # array[int(line[0])][idx] = [int(i) for i in line[1:7]]
                        values = line[1:7]
                        for i in range(6):
                            try:
                                array[i].append(int(values[i]))
                            except KeyError:
                                array[i] = []
                                array[i].append(int(values[i]))


                    except KeyError:
                        # multi
                        # array[int(line[0])] = [0]*len(filenames)
                        # array[int(line[0])][i] = int(line[7])
                        # single
                        # array[int(line[0])] = [0]*7
                        # if line[7] not in line[1:7]:
                        array[int(line[0])] = [int(i)+1 for i in line[1:7]]
                

                    # array[i].append(hexdec)

    # print(array)
    array = np.array(list(array.values()))
    # array = np.rot90(array)
    print(array)
    nmf = NMF(n_components=150, init='random')
    W = nmf.fit_transform(array)
    H = nmf.components_
    # melt = []
    # for key in sorted(array.keys()):
    #     # print(key, strain_dict[key])
    #     for bases in zip(*array[key]):
    #         print(bases)
    #         if sum(bases) != 0:
    #             melt.append(bases)
    # array = np.array(melt)
    # ica = FastICA(n_components=6, max_iter=1000)
    # S_ = ica.fit_transform(array)  # Reconstruct signals
    # A_ = ica.mixing_  # Get estimated mixing matrix
    #
    # pca = PCA(n_components=6)
    # H = pca.fit_transform(array)  # Reconstruct signals based on orthogonal components
    #
    # #############################################################################
    # # Plot results
    #
    # plt.figure()

    models = [array, W, H]
    names = ['Observations (mixed signal)',
             'NMF recovered signals',
             'NMF Clusters']
    colors = ['red', 'steelblue', 'orange', 'pink', 'grey', 'purple']

    for ii, (model, name) in enumerate(zip(models, names), 1):
        plt.subplot(3, 1, ii)
        plt.title(name)
        for sig, color in zip(model, colors):
            plt.plot(sig, 'o', color=color)

    plt.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.46)
    plt.show()
    plt.savefig('strain_unmixed.png')
    # # Mix data


    # pca = PCA()
    # S_pca_ = pca.fit(array).transform(array)
    #
    # ica = FastICA()
    # S_ica_ = ica.fit(array).transform(array)  # Estimate the sources
    #
    # S_ica_ /= S_ica_.std(axis=0)


    # plt.figure()
    # # plt.subplot(2, 2, 1)
    # # plot_samples(S / S.std())
    # # plt.title('True Independent Sources')
    #
    # axis_list = [pca.components_.T, ica.mixing_]
    # plt.subplot(2, 2, 2)
    # plot_samples(array / np.std(array), axis_list=axis_list)
    # legend = plt.legend(['PCA', 'ICA'], loc='upper right')
    # legend.set_zorder(100)
    #
    # plt.title('Observations')
    #
    # plt.subplot(2, 2, 3)
    # plot_samples(S_pca_ / np.std(S_pca_, axis=0))
    # plt.title('PCA recovered signals')
    #
    # plt.subplot(2, 2, 4)
    # plot_samples(S_ica_ / np.std(S_ica_))
    # plt.title('ICA recovered signals')
    #
    # plt.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.36)
    # plt.show()
    return array, W, H


def perform_nmf(filename, k=10):
    # array = [[] for i in filenames]
    array = OrderedDict()
    for (idx, filename) in enumerate(filename):
        with open(filename) as f:
            for line in f:
                if line.startswith('pos'):
                    continue
                else:
                    line = line.strip().split()

                    values = line[1:]
                    array[line[0]] = values

    # print(array)
    col_ids = list(array.keys())
    array = np.array(list(array.values()))
    array = array.T
    lsnmf = nimfa.Lsnmf(array, seed='random_vcol', rank=k, max_iter=10, update='divergence',
                objective='div')
    lsnmf_fit = lsnmf()
    # best_rank = lsnmf.estimate_rank(rank_range=range(8,15))
    # for rank, values in best_rank.items():
    #     print('Rank: %d' % rank)
    #     print('Rss: %5.4f' % values['rss'])
    #     print('Evar: %5.4f' % values['evar'])
    #     print('K-L: %5.4f' % values['kl'])
    # print(best_rank)
    print('Rss: %5.4f' % lsnmf_fit.fit.rss())
    print('Evar: %5.4f' % lsnmf_fit.fit.evar())
    print('K-L divergence: %5.4f' % lsnmf_fit.distance(metric='kl'))
    print('Sparseness, W: %5.4f, H: %5.4f' % lsnmf_fit.fit.sparseness())

    bins = np.array(lsnmf_fit.fit.predict())[0]
    bin_dict = {}
    for (bin_id, contig_name) in zip(bins, col_ids):
        try:
            bin_dict[bin_id].append(contig_name)
        except KeyError:
            bin_dict[bin_id] = [contig_name]

    return lsnmf_fit, bin_dict


def bin_contigs(bin_dict, assembly_file):
    assembly = pyfaidx.Faidx(assembly_file)
    for (bin, contigs) in bin_dict.items():
        with open('bin.'+str(bin)+'.fna', 'w') as f:
            for contig in contigs:
                seq = assembly.fetch(contig, 1, assembly.index[contig].rlen)
                fasta = ">" + seq.name + '\n'
                fasta += seq.seq + '\n'
                f.write(fasta)



def centre(x):
    meanX = np.mean(x, axis=1, keepdims=True)
    return x - meanX, meanX


def covariance(x):
    mean = np.mean(x, axis=1, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean
    return (m.dot(m.T))/n


def whiten(x):
    # covariance matrix
    coVarM = covariance(x)

    # singular value decomposition
    U, S, V = np.linalg.svd(coVarM)

    # calculate diagonal matrix of eigenvalues
    d = np.diag(1.0 / np.sqrt(S))

    # whitening matrix
    whiteM = np.dot(U, np.dot(d, U.T))

    # Project on to whitening matrix
    Xw = np.dot(whiteM, x)

    return Xw, whiteM


def kurtosis(x):
    n = np.shape(x)[0]
    mean = np.sum((x**1)/n) # Calculate the mean
    var = np.sum((x-mean)**2)/n # Calculate the variance
    skew = np.sum((x-mean)**3)/n # Calculate the skewness
    kurt = np.sum((x-mean)**4)/n # Calculate the kurtosis
    kurt = kurt/(var**2)-3
    return kurt, skew, var, mean


def fastIca(signals, alpha=1, thresh=1e-8, iterations=5000):
    m, n = signals.shape
    # Initialize random weights
    W = np.random.rand(m, m)
    for c in range(m):
        w = W[c, :].copy().reshape(m, 1)
        w = w / np.sqrt((w ** 2).sum())
        i = 0
        lim = 100
        while ((lim > thresh) & (i < iterations)):
            # Dot product of weight and signal
            ws = np.dot(w.T, signals)
            # Pass w*s into contrast function g
            wg = np.tanh(ws * alpha).T
            # Pass w*s into g'
            wg_ = (1 - np.square(np.tanh(ws))) * alpha
            # Update weights
            wNew = (signals * wg.T).mean(axis=1) - wg_.mean() * w.squeeze()
            # Decorrelate weights
            wNew = wNew - np.dot(np.dot(wNew, W[:c].T), W[:c])
            wNew = wNew / np.sqrt((wNew ** 2).sum())
            # Calculate limit condition
            lim = np.abs(np.abs((wNew * w).sum()) - 1)

            # Update weights
            w = wNew

            # Update counter
            i += 1
        W[c, :] = w.T
    return W


test, W, H = read_strainm(['tests/test4.tsv'])

bins, ids = perform_nmf(['tests/10_bins_kmer_counts.tsv'], k=30)
bin_contigs(ids, 'tests/10_bins.fna')


Xc, meanX = centre(test)
Xw, whiteM = whiten(Xc)
W = fastIca(Xw, alpha=1, iterations=500)
unMixed = Xw.T.dot(W.T)
unMixed = (unMixed.T - meanX).T

plt.figure()
# plt.plot(unMixed[0], 'ro', unMixed[1], 'bo')
plt.subplot(3, 1, 1)
plt.plot(test[:,-1], 'ro', test[:,-2], 'bo')

plt.subplot(3, 1, 2)
plt.plot(Xw[:,0], 'ro', Xw[:,1], 'bo')

plt.subplot(3, 1, 3)
plt.plot(unMixed[0], 'ro', unMixed[1], 'bo')

plt.show()
plt.savefig('strainm_unmixed.png')

#
# emc2_image = io.imread("tests/crossword.png", as_gray=True)
# ica = FastICA(n_components=817)
# ica.fit(emc2_image)
#
# # reconstruct image with independent components
# emc2_image_ica = ica.fit_transform(emc2_image)
# io.imshow(emc2_image_ica)
# io.show()
# emc2_restored = ica.inverse_transform(emc2_image_ica)
# # io.figure()
# # show image to screen
# io.imshow(emc2_restored)
# io.show()
