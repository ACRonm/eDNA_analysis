import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def edna():
    print('Evaluating eDNA data')

    species_matrix = pd.read_csv('./data/species_matrix.csv')
    species_matrix.head()

    print(species_matrix.head())

    # Plotting the heatmap
    return
