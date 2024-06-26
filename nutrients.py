import pandas as pd
from edna import edna
from edna import plot_nutrients


def nutrients(data):

    print('Evaluating heavy metals data')

    # heavy metals are from cols 22 to 65
    nutrients = data.iloc[:, 15:21]
    site_code = data['Site code']

    nutrients = nutrients.apply(pd.to_numeric, errors='coerce')
    # replace "BDL with imputed value"
    nutrients = nutrients.replace('BDL', 0)

    # add site code first col
    nutrients.insert(0, 'Site code', site_code)

    print(nutrients)

    # average the values from the water and soil concentrations that contain the same site code
    nutrients = nutrients.groupby('Site code').mean()

    genetic_data = pd.read_csv('data/raw_data.csv')

    genetic_data = genetic_data[['Site code', 'Genetic Diversity (shannon)']]

    plot_nutrients(genetic_data, nutrients, "nutrient_concentrations")

    return


def plot_correlation():
    return


if __name__ == "__main__":
    data = pd.read_csv('data/raw_data.csv')
    nutrients(data)
