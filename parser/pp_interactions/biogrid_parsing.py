import pandas as pd
import requests
from tqdm import tqdm
import time
import re
import numpy as np


# import initial data
input_file = 'biogrid.txt'
output_file = 'biogrid_filtered.csv'
final_file = 'biogrid_final.csv'

data = pd.read_csv(input_file, sep="\t")

print("Data read:")
print(data.head())

columns_to_keep = [
    "SWISS-PROT Accessions Interactor A",
    "SWISS-PROT Accessions Interactor B",
    "Score",
    "Qualifications"
]

missing_columns = [col for col in columns_to_keep if col not in data.columns]
if missing_columns:
    print(f"Error: File does not contain these columns: {missing_columns}")
else:
    filtered_data = data[columns_to_keep]

    filtered_data = filtered_data[
        ~(
                filtered_data['SWISS-PROT Accessions Interactor A'].str.contains('-', na=False) |
                filtered_data['SWISS-PROT Accessions Interactor B'].str.contains('-', na=False)
        )
    ]
    filtered_data.reset_index(drop=True, inplace=True)
    # Save required data
    filtered_data.to_csv(output_file, index=False)

    print(filtered_data.info())
    print(f"Data has been saved into {output_file}")

del data


# Parsing names and sequences
data = pd.read_csv(output_file)


def extract_first_id(value):
    if pd.isna(value):
        return None
    return value.split("|")[0].strip()


data["SWISS-PROT Accessions Interactor A"] = data["SWISS-PROT Accessions Interactor A"].apply(extract_first_id)
data["SWISS-PROT Accessions Interactor B"] = data["SWISS-PROT Accessions Interactor B"].apply(extract_first_id)

# Unique IDs
unique_ids = list(set(data["SWISS-PROT Accessions Interactor A"].dropna().tolist() +
                      data["SWISS-PROT Accessions Interactor B"].dropna().tolist()))
print(f"\nlength of unique{len(unique_ids)}")
# Split ID list
num_parts = 5
chunks = [unique_ids[i::num_parts] for i in range(num_parts)]


def get_protein_data_by_id(protein_id):
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}"
    headers = {"Accept": "application/json"}

    try:
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            data = response.json()
            protein_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", "Unknown")
            protein_sequence = data.get("sequence", {}).get("value", "Unknown")
            return protein_name, protein_sequence
        elif response.status_code == 404:
            return "Not Found", "Not Found"
        else:
            return None, None
    except:
        return None, None


for part_num, chunk in enumerate(chunks, start=1):
    protein_names = {}
    protein_sequences = {}

    for protein_id in tqdm(chunk, desc=f"Processing Part {part_num}"):
        protein_name, protein_sequence = get_protein_data_by_id(protein_id)

        if protein_name is not None and protein_sequence is not None:
            protein_names[protein_id] = protein_name
            protein_sequences[protein_id] = protein_sequence

        time.sleep(1)

    output_data = pd.DataFrame({
        "Protein ID": list(protein_names.keys()),
        "Protein Name": list(protein_names.values()),
        "Protein Sequence": list(protein_sequences.values())
    })
    output_data.to_csv(f"p_data_{part_num}.csv", index=False)

print("Processing finished")

intermediate_files = [f"p_data_{i}.csv" for i in range(1, 6)]
protein_datasets = pd.concat([pd.read_csv(file) for file in intermediate_files])


def extract_value(cell):
    try:
        match = re.search(r"value(.*?)}", cell)
        if match:
            raw_value = match.group(1)

            # Clean till first letter or digit
            cleaned_start = re.sub(r"^[^a-zA-Z0-9]*", "", raw_value)

            # Clean after last letter or digit
            cleaned_value = re.sub(r"[^a-zA-Z0-9]*$", "", cleaned_start)

            return cleaned_value
        else:
            return np.nan
    except Exception as e:
        print(f"Error while processing: {cell}. Error: {e}")
        return np.nan


protein_datasets["Protein Name"] = protein_datasets["Protein Name"].apply(extract_value)

protein_name_map = dict(zip(protein_datasets["Protein ID"], protein_datasets["Protein Name"]))
protein_sequence_map = dict(zip(protein_datasets["Protein ID"], protein_datasets["Protein Sequence"]))


def map_name(protein_id):
    return protein_name_map.get(protein_id, "Unknown")


def map_sequence(protein_id):
    return protein_sequence_map.get(protein_id, "Unknown")


data["name_a"] = data["SWISS-PROT Accessions Interactor A"].apply(map_name)
data["sequence_a"] = data["SWISS-PROT Accessions Interactor A"].apply(map_sequence)
data["name_b"] = data["SWISS-PROT Accessions Interactor B"].apply(map_name)
data["sequence_b"] = data["SWISS-PROT Accessions Interactor B"].apply(map_sequence)

data = data.dropna(subset=["SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B"])

# Save final file
data.to_csv(final_file, index=False)

print(f"File saved as {final_file}")
del data

# Formate data
data = pd.read_csv(final_file)
protein_data_a = data[['name_a', 'sequence_a', 'Qualifications']].rename(columns={
    'name_a': 'name',
    'sequence_a': 'content',
    'Qualifications': 'annotation'
})
protein_data_b = data[['name_b', 'sequence_b', 'Qualifications']].rename(columns={
    'name_b': 'name',
    'sequence_b': 'content',
    'Qualifications': 'annotation'
})

# Add required columns
protein_data_a['type'] = 'protein'
protein_data_a['subtype'] = None
protein_data_a['representation_type'] = 'sequence'

protein_data_b['type'] = 'protein'
protein_data_b['subtype'] = None
protein_data_b['representation_type'] = 'sequence'

# Unite data with dropping duplicates
protein_data = pd.concat([protein_data_a, protein_data_b]).drop_duplicates(subset=['name', 'content']).reset_index(drop=True)

# Create interaction_data
interaction_data = data[['name_a', 'name_b', 'Score']].rename(columns={
    'name_a': 'protein_name_a',
    'name_b': 'protein_name_b',
    'Score': 'kd'
})

# Dropping dupclicates
protein_data.drop_duplicates(subset=['name', 'content'], inplace=True)
interaction_data.drop_duplicates(subset=['protein_name_a', 'protein_name_b'], inplace=True)

# Save data
protein_data.to_csv('protein_data_biogrid.csv', index=False)
interaction_data.to_csv('protein_protein_interaction_data.csv', index=False)

print("Protein Data:")
print(protein_data.info())

print("\nInteraction Data:")
print(interaction_data.info())
