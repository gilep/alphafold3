import os
import json
import datetime
from alphafold3 import structure

def update_release_date(mmcif_file: str):
    with open(mmcif_file) as f:
        mmcif_str = f.read()

    DATE = datetime.date(1900, 1, 1)

    struct = structure.from_mmcif(mmcif_str, include_bonds=True)
    struct = struct.copy_and_update_globals(release_date=DATE)
    updated_mmcif = struct.to_mmcif()

    with open(mmcif_file, 'w') as f:
        f.write(updated_mmcif)

def split_cif_by_chains(mmcif_file: str, output_dir: str, chains_list=None):
    """
    Splits an mmCIF file into separate files for each chain and saves them in the output directory.

    :param mmcif_file: Path to the input mmCIF file.
    :param output_dir: Directory where the split chain mmCIF files will be saved.
    :param chains_list: Optional list of chain IDs to save. If None, all chains are saved.
    """
    with open(mmcif_file) as f:
        mmcif_str = f.read()

    struct = structure.from_mmcif(mmcif_str, include_bonds=True)
    chains_struct = struct.split_by_chain()

    os.makedirs(output_dir, exist_ok=True)

    base_filename = os.path.splitext(os.path.basename(mmcif_file))[0]

    for chain_struct in chains_struct:
        chain_id = chain_struct.chains[0]

        # Check if chains_list is provided and if chain_id is in the list
        if chains_list is None or chain_id in chains_list:
            chain_mmcif = chain_struct.to_mmcif()
            output_file = os.path.join(output_dir, f"{base_filename}_chain_{chain_id}.cif")
            with open(output_file, "w") as f:
                f.write(chain_mmcif)

    print(f"Split mmCIF saved in {output_dir}")

def alignment_files_from_json(json_path: str, output_dir: str):
    """
    Extracts unpaired and paired alignment data from an AF3 JSON file and writes them to separate A3M files
    for each sequence based on chain ID (first element if list of IDs).

    :param json_path: Path to the AF3 JSON file.
    :param output_dir: Path to the output directory.
    """
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    base_name = os.path.splitext(os.path.basename(json_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    for entry in af3_data["sequences"]:
        for key, value in entry.items():
            if isinstance(value, dict) and "id" in value:
                chain_id = value["id"]

                if type(chain_id) is list and len(chain_id) > 0:
                    chain_id = chain_id[0]

                if "unpairedMsa" in value:
                    unpaired_file = os.path.join(output_dir, f"{base_name}_unpaired_chain_{chain_id}.a3m")
                    with open(unpaired_file, 'w') as f:
                        f.write(value["unpairedMsa"])

                if "pairedMsa" in value:
                    paired_file = os.path.join(output_dir, f"{base_name}_paired_chain_{chain_id}.a3m")
                    with open(paired_file, 'w') as f:
                        f.write(value["pairedMsa"])


if __name__ == "__main__":
    template_dir = "/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized1"

    for file in os.listdir(template_dir):
        if not file.endswith('.cif') or "chain_" in file:
            continue
        mmcif_path = os.path.join(template_dir, file)
        update_release_date(mmcif_path)
        split_cif_by_chains(mmcif_path, template_dir, chains_list=["A", "B", "C", "D"])
