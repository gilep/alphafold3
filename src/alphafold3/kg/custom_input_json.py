import json
import string

def copy_input_json(json_path: str,  out_path: str, new_name: str,):
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    af3_data['name'] = new_name

    with open(out_path, 'w') as f:
        json.dump(af3_data, f, indent=4)

def add_path_to_msa(json_path: str, chain_id: str, paired_msa_path: str, unpaired_msa_path: str):
    """
    Adds paths for paired and unpaired MSA to a specific chain in the JSON file.

    :param json_path: Path to the AF3 JSON file.
    :param chain_id: Chain ID to update.
    :param paired_msa_path: Path to the paired MSA file.
    :param unpaired_msa_path: Path to the unpaired MSA file.
    """
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    for entry in af3_data["sequences"]:
        for key, value in entry.items():
            if isinstance(value, dict) and "id" in value and value["id"] == chain_id:
                value["pairedMsaPath"] = paired_msa_path
                value["unpairedMsaPath"] = unpaired_msa_path

    with open(json_path, 'w') as f:
        json.dump(af3_data, f, indent=4)

def change_input_json_version(json_path: str, version: int):
    with open(json_path, 'r') as f:
        af3_data = json.load(f)
    af3_data['version'] = version
    with open(json_path, 'w') as f:
        json.dump(af3_data, f, indent=4)

def split_by_chains(json_path: str) :
    """
    Reads a JSON file, splits sequence entries with multiple chain IDs into separate sequences,
    and writes the modified data the same file.

    :param json_path: Path to the AF3 JSON input file.
    """
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    new_sequences = []
    for entry in af3_data["sequences"]:
        for key, value in entry.items():
            if isinstance(value, dict) and "id" in value and isinstance(value["id"], list):
                for chain_id in value["id"]:
                    new_entry = {key: value.copy()}  # Copy original sequence entry
                    new_entry[key]["id"] = chain_id  # Assign individual chain ID
                    new_sequences.append(new_entry)
            else:
                new_sequences.append(entry)  # Keep unchanged entries

    af3_data["sequences"] = new_sequences

    with open(json_path, 'w') as f:
        json.dump(af3_data, f, indent=4)

def glycan_type_to_data(type_glycan):
    if type_glycan is None:
        return None
    elif type_glycan == "NAG":
        glycan_ccds = ['NAG',]
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = None
    elif type_glycan == "NAG(NAG)":
        glycan_ccds = ['NAG', 'NAG']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       ]
    elif type_glycan == "NAG(NAG(BMA))":
        glycan_ccds = ['NAG', 'NAG', 'BMA']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']]
                       ]
    elif type_glycan == "NAG(NAG(BMA(MAN)(MAN)))":
        glycan_ccds = ['NAG', 'NAG', 'BMA', 'MAN', 'MAN']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']],
                       [[3,'O3'],[4, 'C1']],
                       [[3,'O6'],[5, 'C1']],
                       ]
    else:
        raise ValueError(f"Glycan type {type_glycan} is not supported.")
    return glycan_ccds, res_atom_id, glyc_connection_atom, inter_bonds

def next_chain_id(existing_ids):
    """Generate the next available chain ID (A-Z, then AA, AB, etc.)."""
    all_ids = set(existing_ids)

    # Try single-letter IDs first (A-Z)
    for letter in string.ascii_uppercase:
        if letter not in all_ids:
            return letter

    # Extend to two-letter IDs if necessary (AA, AB, AC, ...)
    for first in string.ascii_uppercase:
        for second in string.ascii_uppercase:
            candidate = first + second
            if candidate not in all_ids:
                return candidate
    raise ValueError("Too many chains! No available IDs.")

def generate_indices(query_range: tuple[int]) -> tuple[list[int], list[int]]:
    """
    Generate query and template indices based on the given range.
    :param query_range: Tuple (start, end) for query indices.
    :return: queryIndices and templateIndices lists.
    """
    query_indices = list(range(query_range[0]-1, query_range[1]))
    template_indices = list(range(len(query_indices)))
    return query_indices, template_indices

def add_glycan(json_path, chain_id: str, position: int, type_glycan: str):
    """
    Add glycosylation to the protein chain with the given position.

    :param json_path: Path to the AF3 JSON input file.
    :param type_glycan: Description of the glycan e.g., "NAG" or "NAG(NAG)".
    :param chain_id: Chain ID e.g., "A"
    :param position: Position of the glycosylation. e.g., 42
    """
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    all_chains = []
    for seq in af3_data.get("sequences", []):
        current_chain = list(seq.values())[0]["id"]
        if type(current_chain) is list:
            all_chains.extend(current_chain)
        elif type(current_chain) is str:
            all_chains.append(current_chain)

    glycan_chain_id = next_chain_id(all_chains)
    if glycan_chain_id is None:
        raise ValueError('Chain values are out of range A-Z')

    try:
        # TODO: deal with glycosylation that is not supported
        glycan_ccds, res_atom_id, glyc_connection_atom, inter_bonds = glycan_type_to_data(type_glycan)
    except ValueError as e:
        raise ValueError(e)

    # make new chain
    glycan_chain_data = {"ligand" : {
                            "id" : glycan_chain_id,
                            "ccdCodes" : glycan_ccds
                            }
                        }

    all_chains.append(glycan_chain_id)
    af3_data["sequences"].append(glycan_chain_data)

    # processing bonds
    bonded_atom_pairs = []
    if glyc_connection_atom and res_atom_id:
        bonded_atom_pairs.append([[chain_id, position, res_atom_id],
                                  [glycan_chain_id, *glyc_connection_atom]])

    if inter_bonds:
        bonded_atom_pairs += [[[glycan_chain_id, *atom1], [glycan_chain_id, *atom2]] for atom1, atom2 in inter_bonds]

    if bonded_atom_pairs:
        af3_data["bondedAtomPairs"] = af3_data.get("bondedAtomPairs", []) + bonded_atom_pairs
    else:
        return

    # Update the input file
    with open(json_path, 'w') as f:
        json.dump(af3_data, f, indent=4)
        print(f"Added chain {glycan_chain_id} with glycan {type_glycan} to {json_path}")

def add_protein_template(json_path, chain_id, mmcif_path, query_range):
    """
    Adds a template entry to the protein section of the given AF3 JSON file.

    :param json_path: Path to the AF3 JSON input file.
    :param chain_id: The chain ID of the protein to which the template should be added.
    :param mmcif_path: Path to the mmCIF file.
    :param query_range: Tuple (start, end) for query indices.
    """

    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    for seq in af3_data.get("sequences", []):
        if "protein" in seq and seq["protein"].get("id") == chain_id:
            # Ensure 'templates' key exists
            if "templates" not in seq["protein"]:
                seq["protein"]["templates"] = []


    # Generate indices
    query_indices, template_indices = generate_indices(query_range)

    # Locate the protein section with the specified chain ID
    for seq in af3_data.get("sequences", []):
        if "protein" in seq and seq["protein"].get("id") == chain_id:
            # Ensure 'templates' key exists
            if "templates" not in seq["protein"]:
                seq["protein"]["templates"] = []

            # Add the new template data
            template_data = {
                "mmcifPath": mmcif_path,
                "queryIndices": query_indices,
                "templateIndices": template_indices
            }
            seq["protein"]["templates"].append(template_data)

            # Save the updated JSON file
            with open(json_path, 'w') as f:
                json.dump(af3_data, f, indent=4)

            print(f"Template added to protein chain {chain_id} in {json_path}")
            return

    print(f"Protein chain {chain_id} not found in {json_path}")

# Example usage
if __name__ == "__main__":

    pure_json = "/g/kosinski/kgilep/flu_na_project/na_nc07/af3/input_json/na_nc07.json"

    structure_name = "t2cac4try4"
    json_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/input_json/test/{structure_name}.json"

    copy_input_json(pure_json, json_path, structure_name)

    na_chains = ["A", "B", "C", "D"]
    templates_path_dict = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/templates/try2/T2CAC4_nohead_chain_{id}.mmcif"
                           for id in na_chains}
    paired_msa_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/templates/try2/t2cac4_data_paired_chain_A.a3m"
    unpaired_msa_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/templates/try2/t2cac4_data_unpaired_chain_A.a3m"
    query_range = (1, 63)
    glycosylation_list = [42, 50, 58, 63, 68, 88, 235]
    type_glycan = "NAG(NAG(BMA(MAN)(MAN)))"

    split_by_chains(json_path)
    change_input_json_version(json_path, 2)

    for chain_id in na_chains:
        add_protein_template(json_path, chain_id, templates_path_dict[chain_id], query_range)
        for glycan_num in glycosylation_list:
            add_glycan(json_path, chain_id, glycan_num, type_glycan)
        # add_glycan(json_path, chain_id, 58, 'NAG' )
        add_path_to_msa(json_path, chain_id, paired_msa_path, unpaired_msa_path)

