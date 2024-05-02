from nutrients import nutrients
from heavy_metals import heavy_metals
from edna import edna
import pandas as pd


def menu():
    print("1. Nutrients")
    print("2. Heavy Metals")
    print("3. eDNA")
    return int(input("Enter your choice: "))


def main():
    data = pd.read_csv('data/raw_data.csv', encoding='ISO-8859-1')

    while True:
        try:
            choice = menu()
            if choice == 1:
                nutrients(data)
            elif choice == 2:
                heavy_metals(data)
            elif choice == 3:
                edna(data=None)
            else:
                print("Invalid choice. Please try again.")
                continue
        except Exception as e:
            print("Invalid choice. Please try again.")
            continue


if __name__ == "__main__":
    main()
