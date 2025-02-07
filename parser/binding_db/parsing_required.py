import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import csv

sdf_file = "BindingDB_All_2D.sdf"
output_csv = "binding_db_selected_properties.csv"

required_properties = [
    "BindingDB Ligand Name",
    "Target Name",
    "Target Source Organism According to Curator or DataSource",
    "Ki (nM)",
    "IC50 (nM)",
    "Kd (nM)",
    "EC50 (nM)",
    "pH",
    "Temp C",
    "PubChem CID of Ligand",
    "BindingDB Target Chain Sequence",
    "UniProt (SwissProt) Recommended Name of Target Chain",
    "UniProt (SwissProt) Entry Name of Target Chain",
    "UniProt (SwissProt) Primary ID of Target Chain"
]

# Extracting required properties
suppl = Chem.SDMolSupplier(sdf_file)

with open(output_csv, mode="w", newline="", encoding="utf-8") as csv_file:
    writer = csv.writer(csv_file)

    writer.writerow(required_properties + ["smiles"])

    for mol in tqdm(suppl, desc="Extracting selected properties", unit="molecule"):
        if mol is None:
            continue

        row_data = []
        for prop in required_properties:
            value = mol.GetProp(prop) if mol.HasProp(prop) else None
            row_data.append(value)

        smiles = Chem.MolToSmiles(mol) if mol is not None else None
        row_data.append(smiles)

        writer.writerow(row_data)

print(f"Selected properties have been extracted and saved in {output_csv}")

# Initial data import
df = pd.read_csv('binding_db_selected_properties.csv')
print('data have been imported')

# Create protein data
protein_data = pd.DataFrame({
    "name": df["UniProt (SwissProt) Recommended Name of Target Chain"],
    "type": "protein",
    "subtype": None,
    "representation_type": "sequence",
    "content": df["BindingDB Target Chain Sequence"].str.replace(" ", ""),  # remove spaces
    "annotation": df["Target Source Organism According to Curator or DataSource"]
})
print('protein data have been formed')

# Create molecule data
molecule_data = pd.DataFrame({
    "name": df["BindingDB Ligand Name"],
    "type": "small_molecule",
    "subtype": None,
    "representation_type": "smiles",
    "content": df["smiles"],
    "annotation": None
})
print('molecule data have been formed')

# Create interaction_data
interaction_data = pd.DataFrame({
    "protein_name": df["UniProt (SwissProt) Recommended Name of Target Chain"],
    "small_molecule_name": df["BindingDB Ligand Name"],
    "kd": df["Kd (nM)"]
})
print('interaction data data have been formed')

# Save data
protein_data.to_csv("protein_data.csv", index=False)
molecule_data.to_csv("molecule_data.csv", index=False)
interaction_data.to_csv("interaction_data.csv", index=False)

print("Process finished")
