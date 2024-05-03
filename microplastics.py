import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_microplastics(data):

    print('Evaluating microplastics data')

    print(data.head())

    # barplot for each site
    sns.barplot(x='Site code', y='Number of microplastics',
                data=data, color='#7C9163', edgecolor='black')

    plt.title('Number of microplastics per site')
    plt.xlabel('Site code')
    plt.ylabel('Number of microplastics')
    plt.xticks(rotation=45)

    plt.savefig('./plots/Microplastics/absolute_microplastics.png')
    plt.close()

    # barplot for each site
    sns.barplot(x='Site code', y='Microplastics per litre',
                data=data, color='#7C9163', edgecolor='black')
    plt.title('Microplastics per litre per site')

    plt.xlabel('Site code')
    plt.ylabel('Microplastics per litre')
    plt.xticks(rotation=45)
    # save the plot
    plt.savefig('./plots/Microplastics/relative_microplastics.png')
    plt.close()

    return


def microplastics(data):

    data = data[['Site code', 'Number of microplastics',
                 'Microplastics per litre']]

    # average the values
    data = data.groupby('Site code').mean()

    while True:
        choice = menu()
        if choice == 1:
            plot_microplastics(data)
        elif choice == 2:
            continue
        elif choice == 3:
            break

    return


def menu():

    while (True):
        try:
            print("Choose a microplastics analysis:")

            print("1. Analyse microplastics data")
            print("2. Plot correlation between microplastics and eDNA")
            print("3. Back")
            return int(input("Enter your choice: "))

        except Exception as e:
            print("Invalid choice. Please try again.")
            continue
