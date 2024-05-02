import seaborn as sns
import pandas as pd

df = pd.read_csv('nutrients_vs_ph.csv')

sites = df['Site code']

ph = df['pH']

# nutrients = drop Site code and pH columns
nutrients = df.drop(columns=['Site code', 'pH'])

# for each nutrient col, plot regression against pH (pH on X)
for nutrient in nutrients.columns:
    ph = pd.to_numeric(ph, errors='coerce')
    nutrient_data = pd.to_numeric(nutrients[nutrient], errors='coerce')

    sns.regplot(x=ph, y=nutrient_data)

    # filename = first 3 words
    filename = nutrient.split(' ')[0:3]

    scatterplot = sns.scatterplot(x=ph, y=nutrient_data)

    scatterplot.set_title(f'{nutrient} vs pH')
    scatterplot.set_xlabel('pH')
    scatterplot.set_ylabel(nutrient)

    # set legend
    scatterplot.legend()

    # save the plot as a .png file
    scatterplot.get_figure().savefig(f'{filename}_vs_pH.png')

    # clear
    scatterplot.clear()
