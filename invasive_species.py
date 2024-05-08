import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

species_data = pd.read_csv('data/species_matrix.csv')

water_quality_data = pd.read_csv('data/raw_data.csv')

# select tds, and nutrients cols, filter using regex
water_quality_data = water_quality_data.filter(
    regex='(Site code|TDS ppt|nitrogen|phosphorus|DO|pH|Turbidity)')

# numeric cols are all the cols except the first one
numeric_cols = water_quality_data.columns[1:]

water_quality_data[numeric_cols] = water_quality_data[numeric_cols].apply(
    pd.to_numeric, errors='coerce')


# average all the data with the same Site code
water_quality_data = water_quality_data.groupby(
    'Site code').mean().reset_index()

print(water_quality_data)


carpio = species_data[species_data['Species'] == 'Cyprinus carpio']
carpio_sites = carpio.loc[:, (carpio != 0).any(axis=0)]

# drop the species, and site name columns
carpio_sites = carpio_sites.drop(
    columns=['Species', 'Classification', 'Habitat'])

# drop the first row
carpio_sites = carpio_sites.drop(carpio_sites.index[0])

# carpio sites = the col names
carpio_sites = carpio_sites.columns

# print(carpio_sites)

water_quality_with_carps = water_quality_data[
    water_quality_data['Site code'].isin(carpio_sites)]

water_quality_without_carps = water_quality_data[
    ~water_quality_data['Site code'].isin(carpio_sites)]

water_quality_with_carps = water_quality_with_carps.sort_values(
    'Site code', ascending=True)
water_quality_without_carps = water_quality_without_carps.sort_values(
    'Site code', ascending=True)

# plot the Vis Turbidity
plt.bar(water_quality_with_carps['Site code'],
        water_quality_with_carps['Vis Turbidity'], label='With Carp')
plt.bar(water_quality_without_carps['Site code'],
        water_quality_without_carps['Vis Turbidity'], label='Without Carp')
plt.xlabel('Site Code')
plt.ylabel('Turbidity')
plt.title('Turbidity')
plt.xticks(rotation=45)
plt.legend()
plt.tight_layout()

# save
plt.savefig('turbidity.png')

# mean of the Vis Turbidity
print(water_quality_with_carps['Vis Turbidity'].mean())

print(water_quality_without_carps['Vis Turbidity'].mean())


# create a new dataframe with "Turbidity with carp" and "Turbidity without carp"
turbidity = pd.DataFrame({
    'Turbidity with carp': water_quality_with_carps['Vis Turbidity'],
    'Turbidity without carp': water_quality_without_carps['Vis Turbidity']
})

# transpose

# save the data to a csv file
turbidity.to_csv('turbidity.csv', index=False)

data = pd.read_csv('turbidity.csv')

with_carp = data['Turbidity with carp'].dropna()
without_carp = data['Turbidity without carp'].dropna()

print(with_carp)
print(without_carp)

t_stat, p_val = stats.ttest_ind(with_carp, without_carp)

print(f'The t-statistic is {t_stat} and the p-value is {p_val}')
