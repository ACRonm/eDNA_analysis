import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px


def plot_correlation(corr, data_type):
    # plot the correlation matrix
    annotation = False

    if len(corr.columns) < 10:
        annotation = True

    plt.figure(figsize=(10, 6))
    sns.heatmap(corr, annot=annotation, cmap='coolwarm')
    plt.title('Correlation Matrix')
    plt.xticks(rotation=90)  # Rotate x-axis ticks by 45 degrees
    plt.tight_layout()  # Adjust the layout to prevent tick labels from being cut off
    plt.savefig(f'plots/eDNA/correlation_matrix_{data_type}.png')

    return


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
        plot_species_abundance()

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
        if len(merged_data.columns) > 10:
            # pairplot merged data
            sns.pairplot(data=merged_data,
                         x_vars=merged_data.columns[0:4],
                         y_vars=['Genetic Diversity (shannon)'])
        else:
            for i in range(0, len(merged_data.columns), 4):
                sns.pairplot(data=merged_data,
                             x_vars=merged_data.columns[i:i + 4],
                             y_vars=['Genetic Diversity (shannon)'])
                # save the plot
                plt.savefig(
                    f'plots/eDNA/genetic_diversity_vs_{data_type}_{i}.png')

        # get the correlation between genetic diversity and all other columns
        corr = merged_data.corr()

        # print the corr between genetic diversity and all other columns
        genetic_correlation = corr['Genetic Diversity (shannon)']

        print("Correlation between genetic diversity and all other columns:")

        print(genetic_correlation)

        print("Saving the correlation to a csv file...")
        # save the correlation to a csv file
        genetic_correlation.to_csv(f'data/genetic_correlation_{data_type}.csv')

        plot_correlation(corr, data_type)

        print("Plotting the scatter matrix...")
    else:
        plt.figure(figsize=(10, 6))
        plt.bar(genetic_data['Site code'],
                genetic_data['Genetic Diversity (shannon)'])
        plt.xlabel('Site Code')
        plt.ylabel('Genetic Diversity (Shannon)')
        plt.title('Genetic Diversity by Site')
        plt.xticks(rotation=45)  # Rotate x-axis ticks by 90 degrees
        plt.savefig('./plots/eDNA/genetic_diversity.png')

    return


def plot_species_abundance():

    df = pd.read_csv('data/species_matrix_transposed.csv')

    df.drop(columns=['Site code'], inplace=True)

    df = pd.DataFrame(dict
                      (species=df.columns,
                       abundance=df.sum(axis=0)))

    fig = px.line_polar(df, r='abundance', theta='species', line_close=True)

    fig.update_traces(fill='toself')

    fig.update_layout(polar=dict(radialaxis=dict(type='log')))

    # save the plot
    fig.write_image('plots/eDNA/species_abundance.png', width=1200, height=800)

    fig.show()

    return


if __name__ == '__main__':
    edna(data=None)
