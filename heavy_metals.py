import pandas as pd
import seaborn as sns
import sys
from edna import edna


def heavy_metals(data):

    print('Evaluating heavy metals data')

    # get all columns with heavy metals

    # drop first row
    data = data.drop(0)

    # heavy metals are from cols 22 to 65
    heavy_metals = data.iloc[:, 22:65]
    site_code = data['Site code']

    # add site code first col
    heavy_metals.insert(0, 'Site code', site_code)

    # average the values from the water and soil concentrations that contain the same site code
    heavy_metals = heavy_metals.groupby('Site code').mean()

    print(heavy_metals)

    # water concentrations = cols that contain "Water"
    water_concentrations = heavy_metals.filter(like='Water')
    soil_concentrations = heavy_metals.filter(like='Soil')

    edna(data)

    # average the values from the water and soil concentrations that contain the same site code

    return


def plot_correlation():
    return
