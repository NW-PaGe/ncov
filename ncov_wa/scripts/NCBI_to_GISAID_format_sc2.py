import sys
import re
import os
import pandas as pd
import numpy as np
from Bio import SeqIO

input_metadata = sys.argv[1]
input_sequences = sys.argv[2]

output_metadata = sys.argv[3]
output_sequences = sys.argv[4]

######################
#   Process metadata
######################

#read in metadata
wa_metadata_raw = pd.read_csv(input_metadata, sep="\t", dtype=str)

#create a copy of wa_metadata_raw
wa_metadata_clean = wa_metadata_raw.copy()

### modified this next bit to work as a function:
def clean_date_columns(metadata):
    """
    CLEAN DATE COLUMNS
    args:
        metadata: pd dataframe with dirty columns
    output:
        pd dataframe with clean date columns
    """
    date_raw = wa_metadata_clean["Isolate Collection date"].str.strip().str.replace("/", "-", regex=False)
        #defining this so the operation doens't need to be repeated 6x

    conditions = [
        date_raw.str.match(r"^\d{4}$").fillna(False),           # e.g. "2024"       → "2024-XX-XX"
        date_raw.str.match(r"^\d{4}-\d{2}$").fillna(False),     # e.g. "2024-05"    → "2024-05-XX"
        date_raw.str.match(r"^\d{4}-\d{2}-\d{2}$").fillna(False) # e.g. "2024-05-01" → unchanged
    ]
    choices = [
        date_raw + "-XX-XX",
        date_raw + "-XX",
        date_raw, #if date is already in correct format, no modifications
    ]
    #Formatting 'date' and 'Release_Date' columns
    wa_metadata_clean["date"]         = np.select(conditions, choices, default=None)
    wa_metadata_clean["Release_Date"] = np.select(conditions, choices, default=None)
    return wa_metadata_clean # send only the end result dataframe outside of the function (i.e. 'conditions' disappears)

# now call the function on the data with those filthy date data:
wa_metadata_clean = clean_date_columns(metadata = wa_metadata_clean)

# ---rename columns: original -> new name
wa_metadata_clean = wa_metadata_clean.rename(columns={
    "Accession":                     "strain",
    "Host Name":                     "host",
    "Geographic Region":             "region",
    "SRA Accessions":                "sra_accession",
    "Submitter Names":               "authors",
    "Virus Pangolin Classification": "pangolin_lineage",
})

# ---Exclude non-WA data

    # excluding if 'Geographic Location' has entry like 'USA:Washington, ID/Idaho' or 'USA: WA District of Columbia'
exclude_pattern = r"Washington,?\s*(Idaho|ID|D\.?C\.?|DISTRICT OF COLUMBIA|Bartlesville)|USA:\s*WA,\s*(Dane|DISTRICT)"

    # Include entries that start with "USA: Washington" or "USA: WA", or end with ",WA" / ", WA"
include_pattern = r"^USA:\s*(Washington|WA)(,|$)| WA$|,WA$"

wa_metadata_clean = wa_metadata_clean[
    ~wa_metadata_clean["Geographic Location"].str.contains(exclude_pattern, flags=re.IGNORECASE, regex=True, na=False) & # do not match to exlude_pattern
    wa_metadata_clean["Geographic Location"].str.contains(include_pattern, flags=re.IGNORECASE, regex=True, na=False) # match in include_pattern
]

# ---Parse 'Geographic Location into  'country', 'division', 'location' variables

    # 'country' includes everything before the colon
wa_metadata_clean["country"] = wa_metadata_clean["Geographic Location"].str.extract(r"^(.+?)\s*:")[0].str.strip()

    # 'division' is hard coded to "Washington"
wa_metadata_clean["division"] = "Washington"

    # 'location' is parsed based on format
conditions = [
    # Format: "USA: Washington, <county>" or "USA: Washington,<county>"
    wa_metadata_clean["Geographic Location"].str.match(r"^USA:\s*(Washington|WA),", case=False),
    # Format: "USA: <county>, WA" or "USA: <county>,WA"
    wa_metadata_clean["Geographic Location"].str.match(r"^USA:\s*.+,\s*WA$", case=False),
]
choices = [
    # Extract everything after the first comma
    wa_metadata_clean["Geographic Location"].str.extract(r"^USA:\s*(?:Washington|WA),\s*(.+)$", flags=re.IGNORECASE)[0].str.strip(),
    # Extract the part between ":" and ",WA"
    wa_metadata_clean["Geographic Location"].str.extract(r"^USA:\s*(.+?)\s*,\s*WA$", flags=re.IGNORECASE)[0].str.strip(),
]

wa_metadata_clean["location"] = np.select(conditions, choices, default="?")


# ---create 'virus' variable
wa_metadata_clean["virus"] = np.where(
    wa_metadata_clean["Virus Name"] == "Severe acute respiratory syndrome coronavirus 2",
    "ncov",
    None
)


######################
#   Process fasta
######################

# ---read in fasta file
sequences_raw = SeqIO.to_dict(SeqIO.parse(input_sequences, "fasta"))


# --- clean fasta headers to only ncbi_accession portion

sequences_clean = {}
for header, record in sequences_raw.items():
    match = re.match(r"^(\w+\.\d+)", header)
    clean_name = match.group(1) if match else header
    clean_name = clean_name.rstrip(",") # 'clean_name' is the final cleaned version of the header
    record.id          = clean_name
    record.name        = clean_name
    record.description = clean_name
    sequences_clean[clean_name] = record


# --- make sure metadata and fasta files match
meta_match = wa_metadata_clean[wa_metadata_clean["strain"].isin(sequences_clean.keys())].copy()
seq_match  = [sequences_clean[strain] for strain in meta_match["strain"] if strain in sequences_clean]

# --- output diagnostic messages
print(f"Total sequences in sequences_raw:        {len(sequences_clean)}")
print(f"Total rows in metadata:          {len(wa_metadata_clean)}")
print(f"Metadata rows matching sequences_raw:    {len(meta_match)}")
print(f"Sequences matching metadata:     {len(seq_match)}")
print(f"First few sequences_raw names:           {list(sequences_clean.keys())[:6]}")
print(f"First few metadata strains:      {wa_metadata_clean['strain'].head(6).tolist()}")


# ---write out cleaned metadata

    # Create output directory to write to
os.makedirs(os.path.dirname(output_metadata), exist_ok=True)

meta_match.to_csv(output_metadata, sep="\t", index=False)
print(f"Metadata written to: {output_metadata}")

# --- Write out cleaned fasta

    # Create output directory to write to
os.makedirs(os.path.dirname(output_sequences), exist_ok=True)

SeqIO.write(seq_match, output_sequences, "fasta")
print(f"Sequences written to: {output_sequences}")
