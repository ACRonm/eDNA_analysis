import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def edna(data):

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
        print(data)

        plot_genetic_data(genetic_data, data)

    else:
        print('No data provided. Using default data...')
        plot_genetic_data(genetic_data, data=None)

    # Plotting the heatmap
    return


def plot_genetic_data(genetic_data, data):
    if data is not None:
        # plot the genetic data vs the provided data in a scatter matrix
        print('Plotting genetic data vs provided data...')

        # get the columns of the provided data
        columns = data.columns
        # get the unique values of the site codes
        site_codes = genetic_data['Site code'].unique()
        # create a dictionary to store the data
        data_dict = {}
        # loop through the site codes

        for site_code in site_codes:
            # get the genetic data for the site code
            genetic_datum = genetic_data[genetic_data['Site code'] == site_code]
            # get the data for the site code
            site_data = data[data['Site code'] == site_code]
            # create a dictionary to store the data
            site_data_dict = {}
            # loop through the columns
            for column in columns:
                # get the data for the column
                column_data = site_data[column].values
                # add the data to the dictionary
                site_data_dict[column] = column_data
            # add the data to the dictionary
            data_dict[site_code] = site_data_dict

        # create a dataframe from the dictionary
        data_df = pd.DataFrame(data_dict)
        # create a scatter matrix
        sns.set(style='whitegrid')
        sns.pairplot(data_df, diag_kind='kde')
        plt.show()

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
