import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def edna(data, data_type):

    genetic_data = pd.read_csv('data/raw_data.csv')
    # select only the columns we need

    genetic_data = genetic_data[['Site code', 'Genetic Diversity (shannon)']]

    # drop rows with missing values
    genetic_data = genetic_data.dropna()
    # get all with unique site codes
    genetic_data = genetic_data.drop_duplicates(subset='Site code')

    if data is not None:
        print('Data provided. Using provided data...')
        # read data
        plot_genetic_data(genetic_data, data, data_type)

    else:
        print('No data provided. Using default data...')
        plot_genetic_data(genetic_data, data=None, data_type='eDNA')

    # Plotting the heatmap
    return


def plot_genetic_data(genetic_data, data, data_type):
    if data is not None:
        # plot the genetic data vs the provided data in a scatter matrix
        print('Plotting genetic data vs provided data...')

        print(genetic_data.head())

        print(data.head())

        # impute missing values
        print("Imputing missing data...")
        data = data.fillna(data.mean())

        # plot genetic diversity vs all other columns
        merged_data = pd.merge(genetic_data, data, on='Site code')

        merged_data.drop(columns=['Site code'], inplace=True)

        # plot genetic diversity vs all other columns 4 at a time

        for i in range(0, len(merged_data.columns), 4):
            sns.pairplot(data=merged_data,
                         x_vars=merged_data.columns[i:i + 4],
                         y_vars=['Genetic Diversity (shannon)'])
            # save the plot
            plt.savefig(f'plots/eDNA/genetic_diversity_vs_{data_type}_{i}.png')

        # get the correlation between genetic diversity and all other columns
        corr = merged_data.corr()

        # print the corr between genetic diversity and all other columns
        genetic_correlation = corr['Genetic Diversity (shannon)']

        print("Correlation between genetic diversity and all other columns:")

        print(genetic_correlation)

        print("Saving the correlation to a csv file...")
        # save the correlation to a csv file
        genetic_correlation.to_csv(f'data/genetic_correlation_{data_type}.csv')

        print("Plotting the scatter matrix...")
    else:

        plt.figure(figsize=(10, 6))
        plt.bar(genetic_data['Site code'],
                genetic_data['Genetic Diversity (shannon)'])
        plt.xlabel('Site Code')
        plt.ylabel('Genetic Diversity (Shannon)')
        plt.title('Genetic Diversity by Site')
        plt.xticks(rotation=45)  # Rotate x-axis ticks by 90 degrees
        plt.savefig('plots/eDNA/genetic_diversity.png')

        plt.show()
        # save plots/eDNA/genetic_diversity.png
    return


if __name__ == '__main__':
    edna(data=None)
