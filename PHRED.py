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
                 color=plt.cm.tab20b(np.arange(len(df.columns))), edgecolor='black')

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

    # plot a pie chart for all the samples
    ax = df.sum().plot(kind='pie', figsize=(10, 10), autopct='%1.1f%%',
                       colors=plt.cm.tab20b(np.arange(len(df.columns))), startangle=90, wedgeprops=dict(width=0.5), textprops=dict(color='black', fontsize=10))
    ax.set_title('Taxonomic Abundance')
    plt.tight_layout()

    # Adjust the position of the labels
    ax.legend(bbox_to_anchor=(1, 0.5), loc='center left')

    # Offset the percentages to avoid overlap
    plt.subplots_adjust(left=0.1, right=0.7)

    plt.savefig('./plots/taxonomic_abundance_pie.png', bbox_inches='tight')

    # calculate the mean Alpha diversity using the Shannon index


def plot_abundance_per_species():
    file_name = './data/species_matrix_transposed.csv'
    df = pd.read_csv(file_name, index_col=0)
    df.index = ['ETT.' + str(i).zfill(2) for i in range(1, len(df.index) + 1)]

    # Remove any rows that have all zeros
    df = df.loc[(df != 0).any(axis=1)]

    # Calculate the total number of reads per sample
    total_reads = df.sum(axis=1)

    # Calculate the relative abundance of each species
    relative_abundance = df.div(total_reads, axis=0)

    # Calculate the mean relative abundance of each species
    mean_relative_abundance = relative_abundance.mean()

    # Sort the mean relative abundance in descending order
    mean_relative_abundance = mean_relative_abundance.sort_values(
        ascending=False)

    # Plot the mean relative abundance of each species
    ax = mean_relative_abundance.plot(kind='bar', figsize=(
        10, 5), color='skyblue', edgecolor='black')

    plt.xlabel('Species')
    plt.ylabel('Mean Relative Abundance')
    plt.title('Mean Relative Abundance of Each Species')

    plt.xticks(rotation=45, ha='right')

    # Add data labels to the bars
    for i, v in enumerate(mean_relative_abundance):
        ax.text(i, v, str(round(v, 2)), ha='center', va='bottom')

    plt.savefig('./plots/mean_relative_abundance.png', bbox_inches='tight')
    plt.close()


def plot_alpha_diversity():
    return


if __name__ == '__main__':

    # plot_reads_per_sample()
    # plot_alpha_diversity()
    plot_abundance_per_species()
    directory = './data/ETT_filtered_sequences/'

    total_qualities = []

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

    mean_quality = np.mean(total_qualities)
    print('Mean quality: ', mean_quality)
