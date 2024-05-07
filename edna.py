import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import statsmodels.api as sm
from statsmodels.formula.api import ols


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

    # average duplicate site codes
    genetic_data = genetic_data.groupby('Site code').mean().reset_index()

    # drop ETT24.12
    genetic_data = genetic_data[genetic_data['Site code'] != 'ETT24.12']

    if data is not None:
        print('Data provided. Using provided data...')
        # read data
        plot_genetic_data(genetic_data, data, data_type)

        if "metal" in data_type:
            plot_metal_relationship(genetic_data, data, data_type)

        plot_nutrients*(genetic_data, data, data_type)

    else:
        # print('No data provided. Using default data...')
        # plot_genetic_data(genetic_data, data=None, data_type='eDNA')
        # plot_species_abundance()
        # plot_stacked_bar_abundance()

        plot_nutrients(genetic_data, data=None, data_type='Nutrients')

    # Plotting the heatmap
    return


def plot_genetic_data(genetic_data, data, data_type):
    if data is not None:
        # plot the genetic data vs the provided data in a scatter matrix
        print('Plotting genetic data vs provided data...')

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

        print("Saving the correlation to a csv file...")

        # sort by correlation
        genetic_correlation = genetic_correlation.sort_values(ascending=False)

        # save the correlation to a csv file
        genetic_correlation.to_csv(f'data/genetic_correlation_{data_type}.csv')

        plot_correlation(corr, data_type)

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

    return


def plot_stacked_bar_abundance():
    df = pd.read_csv('data/species_matrix_transposed.csv')

    # get totals

    # for each value in the site code column, replace "Site_" with "ETT24."
    df['Site code'] = df['Site code'].replace("Site_", "ETT24.", regex=True)

    # plot the stacked bar chart with the species abundance on the y-axis and the site code on the x-axis
    df.set_index('Site code', inplace=True)

    ax = df.plot(kind='bar', stacked=True, figsize=(
        10, 6), colormap='tab20b', edgecolor='black')
    plt.title('Species abundance by site code')
    plt.xlabel('Species')
    plt.ylabel('Abundance')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    # Add bbox_inches='tight' to include legend
    plt.savefig('plots/eDNA/stacked_bar_abundance.png', bbox_inches='tight')

    # get transpose of the dataframe
    df = df.transpose()

    # get totals
    df['Total'] = df.sum(axis=1)

    # sort the dataframe by the total column
    df = df.sort_values(by='Total', ascending=False)

    print(df['Total'])

    return


def plot_metal_relationship(genetic_data, data, data_type):

    # close all plots
    plt.close('all')

    # plot the relationship between the genetic diversity and the metal concentrations
    # merge the genetic data with the metal data
    data = data.filter(regex='Pb|As|Ca|Zn|Ni|Cu|Li|U')
    merged_data = pd.merge(genetic_data, data, on='Site code')

    # drop na
    merged_data.dropna(inplace=True)

    # drop the site code column
    merged_data.drop(columns=['Site code'], inplace=True)

    # perform two-way anova between genetic diversity and metal concentrations
    print("Performing two-way ANOVA between genetic diversity and metal concentrations...")
    # create a formula
    # Assuming merged_data is your DataFrame and has been correctly prepared

    # Create a formula by wrapping column names with spaces in Q("")
    formula = 'Q("Genetic Diversity (shannon)") ~ ' + \
        ' + '.join([f'Q("{col}")' for col in merged_data.columns[1:]])

    print(formula)

    # perform the two-way anova
    model = ols(formula, data=merged_data).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    print(anova_table)

    # save the anova table to a csv file
    anova_table.to_csv(f'data/anova_table_{data_type}.csv')

    # plot linear regression between genetic diversity and metal concentrations
    for i in range(1, len(merged_data.columns)):
        sns.lmplot(x=merged_data.columns[i], y='Genetic Diversity (shannon)',
                   data=merged_data, height=6, aspect=1.5, line_kws={'color': 'red'}, ci=None)

        plt.savefig(
            f'plots/eDNA/genetic_diversity_vs_{data_type}_{i}.png')

    return


def plot_nutrients(genetic_data, data, data_type):

    # multi variable regression
    # merge the genetic data with the nutrient data
    data = data.filter(regex='phosphorus|nitrogen')

    merged_data = pd.merge(genetic_data, data, on='Site code')
    # drop na
    merged_data.dropna(inplace=True)

    # Get unique 'Site code' values
    site_codes = merged_data['Site code'].unique()

    for i in range(1, len(merged_data.columns)):
        plt.figure(figsize=(9, 6))  # Adjust size as needed
        for site in site_codes:
            # Subset the data for this site
            subset = merged_data[merged_data['Site code'] == site]
            # Make a scatter plot with custom color
            sns.scatterplot(x=subset[merged_data.columns[i]],
                            y=subset['Genetic Diversity (shannon)'], label=site)
        # Add a regression line for all data points
        sns.regplot(x=merged_data[merged_data.columns[i]],
                    y=merged_data['Genetic Diversity (shannon)'], scatter=False, color='red')
        plt.xlabel(merged_data.columns[i])
        plt.ylabel('Genetic Diversity (shannon)')
        plt.title(f'Genetic Diversity vs {merged_data.columns[i]}')
        plt.legend(title='Site code')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f'plots/eDNA/genetic_diversity_vs_{data_type}_{i}.png')
        plt.close()

    # drop the site code column
    merged_data.drop(columns=['Site code'], inplace=True)

    # perform two-way anova between genetic diversity and nutrient concentrations
    print("Performing two-way ANOVA between genetic diversity and nutrient concentrations...")
    # create a formula
    # Assuming merged_data is your DataFrame and has been correctly prepared

    # Create a formula by wrapping column names with spaces in Q("")
    formula = 'Q("Genetic Diversity (shannon)") ~ ' + \
        ' + '.join([f'Q("{col}")' for col in merged_data.columns[1:]])

    print(formula)

    # perform the two-way anova
    model = ols(formula, data=merged_data).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    print(anova_table)

    # save the anova table to a csv file
    anova_table.to_csv(f'data/anova_table_{data_type}.csv')

    return


if __name__ == '__main__':
    edna(data=None)
