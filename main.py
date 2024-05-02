from Nutrients_analysis.nutrients import nutrients
from Heavy_metals_analysis.heavy_metals import heavy_metals
from eDNA_analysis.edna import edna


def menu():
    print("1. Nutrients")
    print("2. Heavy Metals")
    print("3. eDNA")
    return int(input("Enter your choice: "))


def main():
    while True:
        choice = menu()
        if choice == 1:
            nutrients()
        elif choice == 2:
            heavy_metals()
        elif choice == 3:
            edna()
        else:
            print("Invalid choice")


if __name__ == "__main__":
    main()
