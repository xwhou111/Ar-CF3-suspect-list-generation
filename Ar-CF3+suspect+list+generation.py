# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
import csv
import itertools
from collections import Counter

SCAFFOLDS = [
    ("benzene", {"C": 6}, 6),
]

SUBSTITUENTS = {
    "CH3": {"C": 1, "H": 3},
    "C2H5": {"C": 2, "H": 5},
    "CHO": {"C": 1, "H": 1, "O": 1},
    "CN": {"C": 1, "N": 1},
    "COOH": {"C": 1, "O": 2, "H": 1},
    "NH2": {"N": 1, "H": 2},
    "NHCH3": {"N": 1, "C": 1, "H": 4},
    "N(CH3)2": {"N": 1, "C": 2, "H": 6},
    "NO2": {"N": 1, "O": 2},
    "OH": {"O": 1, "H": 1},
    "SH": {"S": 1, "H": 1},
    "SCH3": {"S": 1, "C": 1, "H": 3},
    "SO2H": {"S": 1, "O": 2, "H": 1},
    "SO3H": {"S": 1, "O": 3, "H": 1},
    "PO3H2": {"P": 1, "O": 3, "H": 2},
    "F": {"F": 1},
    "Cl": {"Cl": 1},
    "Br": {"Br": 1},
}

MASSES = {
    "C": 12.00000,
    "H": 1.00783,
    "N": 14.00307,
    "O": 15.99491,
    "F": 18.99840,
    "P": 30.97376,
    "S": 31.97207,
    "Cl": 34.96885,
    "Br": 78.91834,
}

ELECTRON_MASS = 0.000549
POLARITIES = {"positive": 1, "negative": -1}
MIN_MZ, MAX_MZ = 0, 400

VARIANTS = {
    "base_ion": {},
    "H_loss": {"H": -1},
    "H2_loss": {"H": -2},
    "HF_loss": {"H": -1, "F": -1},
}

OUTPUT_FILE = "Ar_CF3_suspect_list.csv"


def formula_text(atoms):
    order = ["C", "H"] + sorted(set(atoms) - {"C", "H"})
    return "".join(
        element + (str(atoms[element]) if atoms[element] != 1 else "")
        for element in order
        if atoms.get(element, 0) > 0
    )


def generate():
    results = {}

    for scaffold, core, original_h in SCAFFOLDS:
        for cf3_count in (1, 2):
            for other_count in range(original_h - cf3_count + 1):
                for groups in itertools.combinations_with_replacement(
                    SUBSTITUENTS, other_count
                ):
                    parent = Counter(core)
                    parent.update({
                        "C": cf3_count,
                        "F": 3 * cf3_count,
                        "H": original_h - cf3_count - other_count,
                    })

                    for group in groups:
                        parent.update(SUBSTITUENTS[group])

                    group_text = "+".join(groups) or "none"

                    for polarity, charge in POLARITIES.items():
                        for variant, changes in VARIANTS.items():
                            ion = parent.copy()
                            ion["H"] += charge
                            ion.update(changes)

                            if any(n < 0 for n in ion.values()):
                                continue

                            mz = sum(
                                MASSES[element] * n
                                for element, n in ion.items()
                            ) - charge * ELECTRON_MASS

                            if not MIN_MZ <= mz <= MAX_MZ:
                                continue

                            formula = formula_text(ion)
                            key = (formula, charge)
                            source = (
                                f"{scaffold}; CF3={cf3_count}; "
                                f"groups={group_text}; variant={variant}"
                            )

                            if key not in results:
                                results[key] = {
                                    "ion_formula": formula,
                                    "polarity": polarity,
                                    "charge": charge,
                                    "mz": mz,
                                    "sources": set(),
                                }

                            results[key]["sources"].add(source)

    return sorted(
        results.values(),
        key=lambda row: (row["polarity"], row["mz"]),
    )


def main():
    results = generate()

    with open(
        OUTPUT_FILE, "x", newline="", encoding="utf-8-sig"
    ) as file:
        writer = csv.writer(file)
        writer.writerow([
            "ion_formula", "polarity", "charge", "mz", "sources"
        ])

        for row in results:
            writer.writerow([
                row["ion_formula"],
                row["polarity"],
                row["charge"],
                f'{row["mz"]:.5f}',
                " | ".join(sorted(row["sources"])),
            ])

    print(f"Generated {len(results)} unique candidate ions.")
    print(f"Saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()