#####################
# input json file functions
#####################

import pandas as pd
import json
import os
from typing import List, Dict
import pickle
import gzip

def prepare_crosslink_pickle(respair_path: str, output_path: str, combo_list: List[str], FDR_VALUE: float = 0.02):
    """
    Prepare and save a crosslink dictionary as a compressed pickle file.

    Args:
        respair_path (str): Path to the CSV file containing crosslinking data.
        output_path (str): Path where the pickle file will be saved.
        combo_list (List[str]): List of protein combinations.
    """
    # Load CSV efficiently
    respair_df = pd.read_csv(respair_path, dtype=str).fillna('')

    # Precompute mappings for quick lookup
    alpha_prot_map = respair_df['Alpha protein mapping(s)'].str.split("; ").tolist()
    beta_prot_map = respair_df['Beta protein mapping(s)'].str.split("; ").tolist()
    alpha_pos_map = respair_df['Alpha protein(s) position(s)'].str.split(";").tolist()
    beta_pos_map = respair_df['Beta protein(s) position(s)'].str.split(";").tolist()

    cross_full_dict = {}

    for combo in combo_list:
        prot1, prot2 = combo.split('_')
        cross_set = set()

        for row_idx in range(len(respair_df)):
            match1, match2, pos1, pos2 = False, False, -1, -1

            # Alpha protein matching
            for num1, row_prot in enumerate(alpha_prot_map[row_idx]):
                prot_names = row_prot.split('|')[1:]
                if prot1 in prot_names:
                    match1 = True
                    pos1 = int(alpha_pos_map[row_idx][num1])
                elif prot2 in prot_names:
                    match2 = True
                    pos2 = int(alpha_pos_map[row_idx][num1])

            if not (match1 or match2):
                continue

            # Beta protein matching
            for num2, row_prot in enumerate(beta_prot_map[row_idx]):
                prot_names = row_prot.split('|')[1:]
                if match2 and prot1 in prot_names:
                    match1 = True
                    pos1 = int(beta_pos_map[row_idx][num2])
                if match1 and prot2 in prot_names:
                    match2 = True
                    pos2 = int(beta_pos_map[row_idx][num2])

            if match1 and match2:
                cross_set.add((pos1, pos2))

        if cross_set:
            cross_list = [(pos1 - 1, pos2 - 1, FDR_VALUE) for pos1, pos2 in cross_set]
            cross_full_dict.setdefault(prot1, {})[prot2] = cross_list
        else:
            print(f'{combo} is empty')

    # Save as compressed pickle
    with gzip.open(output_path, 'wb') as f:
        pickle.dump(cross_full_dict, f)

def parse_fasta(fasta_path: str) -> dict:
    """
    Parse a FASTA file to extract protein ID and sequence.

    Args:
        fasta_path (str): Path to the FASTA file.

    Returns:
        dict: A dictionary containing 'sequence'.
    """
    with open(fasta_path, 'r', encoding='utf-8') as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]
        if not lines or lines[0][0] != '>':
            raise ValueError(f"Invalid FASTA format in file: {fasta_path}")
        sequence = ''.join(lines[1:])
    return {"sequence": sequence}

def pkl_to_json_cross(
    combo: List[str],
    cross_dict: Dict[str, Dict[str, List[List[int]]]],
    crosslinker: str = "DSSO"
) -> dict:
    residue_pairs = []
    protein_ids = {protein: chr(65 + idx) for idx, protein in enumerate(combo)}

    for protein1 in combo:
        for protein2 in combo:
            if protein1 == protein2:
                continue
            for direction in [(protein1, protein2), (protein2, protein1)]:
                cross_set = cross_dict.get(direction[0], {}).get(direction[1], [])
                for cross in cross_set:
                    if len(cross) < 2:
                        continue  # Skip invalid entries
                    pos1, pos2 = cross[:2]  # Ignore extra fields if present
                    cross_note = [[protein_ids[direction[0]], pos1 + 1],
                                  [protein_ids[direction[1]], pos2 + 1]]
                    if cross_note not in residue_pairs:
                        residue_pairs.append(cross_note)

    return {
        "name": crosslinker,
        "residue_pairs": residue_pairs
    }

def prepare_alphafold3_json(fasta_paths: List[str], cross_dict: dict, model_seeds: List[int] = [1, 2, 3, 4], version: int = 1) -> dict:
    """
    Prepare a JSON object for AlphaFold3 input to predict heterodimers.

    Args:
        fasta_paths (List[str]): List of paths to protein FASTA files.
        model_seeds (List[int]): List of model seeds for the prediction.

    Returns:
        dict: The prepared JSON object.
    """
    sequences = []
    protein_ids = []

    # Parse each FASTA file and add to sequences list
    for index, fasta_path in enumerate(fasta_paths):
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        protein = parse_fasta(fasta_path)
        protein_id = chr(65 + index)  # Assign IDs starting from 'A'
        protein["id"] = protein_id
        protein_ids.append(os.path.basename(fasta_path).split('.')[0])
        sequences.append({"protein": protein})

    # Construct the job name
    job_name = '_'.join(protein_ids)

    # Add crosslinks
    crosslinks = [pkl_to_json_cross(protein_ids, cross_dict)]

    # Construct the JSON structure
    input_json = {
        "name": job_name,
        "modelSeeds": model_seeds,
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": version,
        "crosslinks": crosslinks
    }
    return input_json

if __name__ == "__main__":
    respair_path = "/g/kosinski/kgilep/flu_hosthost_project/crosslinks/ResPair_2FDR_Decoy_allEnr.csv"

    combo_file = "/g/kosinski/kgilep/flu_hostpat_cross/sample_sheet.txt"
    with open(combo_file, 'r') as f:
       combo_list = [line.strip() for line in f.readlines()]

    crosslinks_full_path = '/g/kosinski/kgilep/flu_hostpat_cross/alphalink/input/2crosslinks.pkl.gz'

    prepare_crosslink_pickle(respair_path, crosslinks_full_path, combo_list)

    with gzip.open(crosslinks_full_path, 'rb') as f:
        cross_dict = pickle.load(f)

    fasta_dir = '/scratch/kgilep/flu_hostpat_cross/af3/data'
    fasta1 = os.path.join(fasta_dir, 'M1.fasta')
    fasta2 = os.path.join(fasta_dir, 'P60709.fasta')
    result = prepare_alphafold3_json([fasta1, fasta2], cross_dict)

    print(json.dumps(result, indent=2))
