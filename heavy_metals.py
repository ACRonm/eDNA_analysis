import pandas as pd
from edna import edna
import matplotlib.pyplot as plt
import numpy as np


def heavy_metals(data):

    print('Evaluating heavy metals data')

    # heavy metals are from cols 22 to 65
    heavy_metals = data.iloc[:, 21:65]
    site_code = data['Site code']

    # add site code first col
    heavy_metals.insert(0, 'Site code', site_code)

    # average the values from the water and soil concentrations that contain the same site code
    heavy_metals = heavy_metals.groupby('Site code').mean()

    # water concentrations = cols that contain "Water"
    water_concentrations = heavy_metals.filter(like='Water')
    soil_concentrations = heavy_metals.filter(like='Soil')

    # plot the heavy metal concentrations
    plot_heavy_metal_concentrations(
        water_concentrations, "water_metal_concentrations")
    plot_heavy_metal_concentrations(
        soil_concentrations, "soil_metal_concentrations")

    print("Analysing correlation between water concentrations and eDNA")
    edna(water_concentrations, "water_metal_concentrations")

    print("Analysing correlation between soil concentrations and eDNA")
    edna(soil_concentrations, "soil_metal_concentrations")

    return


def plot_correlation():
    return


def plot_heavy_metal_concentrations(metal_data, datatype):

    # close any existing plots
    plt.close()

    # remove any rows with NaN values
    metal_data.dropna(inplace=True)

    # filter mg, as, cu, ni, pb, zn with regex
    metal_data = metal_data.filter(regex='(Ca)')

    print(metal_data.head())

    # bar chart
    metal_data.plot(kind='bar', figsize=(20, 10))

    plt.xlabel('Site code')
    plt.ylabel('Concentration (mg/kg)')
    plt.title('Heavy metal concentrations')
    plt.tight_layout()

    # save the plot
    plt.savefig('./plots/Heavy_metals/' + datatype + '.png')
    return


if __name__ == "__main__":
    data = pd.read_csv('data/raw_data.csv')
    heavy_metals(data)
