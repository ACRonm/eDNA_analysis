import pandas as pd
from edna import edna


def heavy_metals(data):

    print('Evaluating heavy metals data')

    # heavy metals are from cols 22 to 65
    heavy_metals = data.iloc[:, 21:65]
    site_code = data['Site code']

    print(heavy_metals.describe())

    # add site code first col
    heavy_metals.insert(0, 'Site code', site_code)

    # average the values from the water and soil concentrations that contain the same site code
    heavy_metals = heavy_metals.groupby('Site code').mean()

    # water concentrations = cols that contain "Water"
    water_concentrations = heavy_metals.filter(like='Water')
    soil_concentrations = heavy_metals.filter(like='Soil')

    print("Analysing correlation between water concentrations and eDNA")
    edna(water_concentrations, "water_metal_concentrations")

    print("Analysing correlation between soil concentrations and eDNA")
    edna(soil_concentrations, "soil_metal_concentrations")

    return


def plot_correlation():
    return


if __name__ == "__main__":
    data = pd.read_csv('data/raw_data.csv')
    heavy_metals(data)
