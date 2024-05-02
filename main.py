from Nutrients_analysis.nutrients import nutrients
from Heavy_metals_analysis.heavy_metals import heavy_metals
from eDNA_analysis.edna import edna
import pandas as pd


def menu():
    print("1. Nutrients")
    print("2. Heavy Metals")
    print("3. eDNA")
    return int(input("Enter your choice: "))


def main():
    data = pd.read_csv('data/raw_data.csv', encoding='ISO-8859-1')

    while True:
        choice = menu()
        if choice == 1:
            nutrients(data)
        elif choice == 2:
            heavy_metals(data)
        elif choice == 3:
            edna(data)
        else:
            print("Invalid choice")


if __name__ == "__main__":
    main()
