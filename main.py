from nutrients import nutrients
from heavy_metals import heavy_metals
from edna import edna
from microplastics import microplastics
import pandas as pd


def menu():
    print("1. Nutrients")
    print("2. Heavy Metals")
    print("3. eDNA")
    print("4. Microplastics")
    print("5. Exit")
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
                edna(data=None, data_type=None)
            elif (choice == 4):
                microplastics(data)
            elif (choice == 5):
                break
            else:
                print("Invalid choice. Please try again.")
                continue
        except Exception as e:
            print(e)
            print("Invalid choice. Please try again.")
            continue


if __name__ == "__main__":
    main()
