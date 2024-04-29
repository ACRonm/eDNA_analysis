import re
import os
import math
import pylab as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
from Bio import SeqIO
import gzip


def plot_fastq_qualities(filename, ax=None, limit=1000):

    fast_parser = SeqIO.parse(gzip.open(filename, 'rt'), 'fastq')

    res = []
    c = 0

    for record in fast_parser:
        c += 1
        res.append(record.letter_annotations['phred_quality'])
        if c > limit:
            break

    df = pd.DataFrame(res)
    l = len(df.columns)

    if ax == None:
        f, ax = plt.subplots(figsize=(7, 5))
    rect = patches.Rectangle(
        (0, 0), l, 20, linewidth=0, facecolor='r', alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0, 20), l, 8, linewidth=0,
                             facecolor='yellow', alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle(
        (0, 28), l, 12, linewidth=0, facecolor='g', alpha=.4)
    ax.add_patch(rect)
    df.mean().plot(ax=ax, c='black')
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    df.plot(kind='box', ax=ax, grid=False, showfliers=False,
            color=dict(boxes='black', whiskers='black'))
    ax.set_xticks(np.arange(0, l, 10))
    ax.set_xticklabels(np.arange(0, l, 10), rotation=45)
    ax.set_xlabel('position(bp)')
    ax.set_xlim((0, l))
    ax.set_ylim((0, 40))

    mean_quality = df.mean().mean()

    print('Mean quality: ', mean_quality
          )
    plt.tight_layout()  # Crop white space around the plot
    return mean_quality


def plot_reads_per_sample():
    filename = './data/species_matrix_transposed.csv'

    # Set the first column as the index
    df = pd.read_csv(filename, index_col=0)
    df.index = ['ETT.' + str(i).zfill(2) for i in range(1, len(df.index) + 1)]

    print(df)

    # Remove any rows that have all zeros
    df = df.loc[(df != 0).any(axis=1)]

    print(df)

    print(len(df.columns))

    ax = df.plot(kind='bar', stacked=True, figsize=(10, 5),
                 color=plt.cm.viridis(np.arange(len(df.columns))), edgecolor='black')

    plt.yscale('linear')
    plt.xlabel('Sampling Site')
    plt.ylabel('Number of Reads per Sample (log)')

    # rotate xticks parallel to the axis
    # Adjust the rotation and alignment of xticks
    plt.xticks(rotation=45, ha='right')

    # move the legend outside the chart
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # save
    plt.savefig('./plots/taxonomic_abundance.png', bbox_inches='tight')
    plt.close()

    # calculate the mean Alpha diversity using the Shannon index


if __name__ == '__main__':

    plot_reads_per_sample()
    directory = './data/ETT_Filtered_sequences/'

    # total_qualities = []

    # if not os.path.exists('./plots'):
    #     os.makedirs('./plots')

    # for filename in os.listdir(directory):
    #     print(filename)
    #     if filename.endswith('.fastq.gz'):
    #         filename = os.path.join(directory, filename)
    #         mean = plot_fastq_qualities(filename)

    #         # save the plot
    #         plt.savefig(
    #             './plots/' + os.path.basename(filename).replace('.fastq.gz', '_quality.png'))

    #         plt.close()

    #         # add mean to list
    #         total_qualities.append(mean)

    # mean_quality = np.mean(total_qualities)
    # print('Mean quality: ', mean_quality)
