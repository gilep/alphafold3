import json
import string
import os

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

def clear_templates_for_chain(chain_id, input_json_path, output_json_path):
    """
    Sets the 'templates' field to null for a specific protein chain in an AlphaFold3 input JSON.

    Parameters:
    - chain_id: the 'id' of the protein to update (e.g., 'A')
    - input_json_path: path to the input JSON file.
    - output_json_path: path to save the modified JSON.
    """
    with open(input_json_path, 'r') as f:
        data = json.load(f)

    for entry in data.get('sequences', []):
        protein = entry.get('protein')
        if protein and protein.get('id') == chain_id:
            protein['templates'] = []

    with open(output_json_path, 'w') as f:
        json.dump(data, f, indent=2)

def change_seeds(json_path: str, num_seeds: int):
    with open(json_path, 'r') as f:
        af3_data = json.load(f)
    af3_data['modelSeeds'] = list(range(1, num_seeds+1))
    with open(json_path, 'w') as f:
        json.dump(af3_data, f, indent=4)

def trim_a3m_lines(position_range, a3m_lines):
    """
    Trims A3M lines based on alignment positions (1-based), counting only uppercase and '-' characters.

    Parameters:
    - position_range: (start, end) — 1-based inclusive positions in the aligned region.
    - a3m_lines: list of lines from an A3M file (with headers and sequences).

    Returns:
    - list of trimmed A3M lines.
    """
    start_pos, end_pos = position_range
    start_pos -= 1
    output_lines = []

    for i, line in enumerate(a3m_lines):
        if line.startswith('>'):
            output_lines.append(line)
        else:
            seq = line.strip()
            pos = 0
            for idx, c in enumerate(seq):
                if c.isupper() or c == '-':
                    pos += 1
                if start_pos == pos:
                    start_trim = idx
                elif end_pos == pos:
                    end_trim = idx
            trimmed_seq = seq[start_trim:end_trim]
            output_lines.append(trimmed_seq + '\n')

    return output_lines

def trim_a3m_file(position_range, input_path, output_path):
    with open(input_path, 'r') as f:
        lines = f.readlines()

    trimmed_lines = trim_a3m_lines(position_range, lines)

    with open(output_path, 'w') as f:
        f.writelines(trimmed_lines)

def trim_protein_chain(protein_range, chain_id, input_json_path, output_json_path):
    start, end = protein_range

    def trim_a3m_string(a3m_string, position_range):
        """Trim A3M block from a string using the same logic as A3M files."""
        lines = a3m_string.strip().splitlines(keepends=False)
        trimmed_lines = trim_a3m_lines(position_range, [line + '\n' for line in lines])
        return ''.join(trimmed_lines)

    with open(input_json_path, 'r') as f:
        data = json.load(f)

    for entry in data['sequences']:
        protein = entry.get('protein')
        if not protein or protein.get('id') != chain_id:
            continue

        # Trim sequence
        full_seq = protein['sequence']
        protein['sequence'] = full_seq[start - 1:end]

        # Trim modifications
        if 'modifications' in protein:
            trimmed_mods = []
            for mod in protein['modifications']:
                pos = mod['ptmPosition']
                if start <= pos <= end:
                    mod_copy = mod.copy()
                    mod_copy['ptmPosition'] = pos - start + 1
                    trimmed_mods.append(mod_copy)
            protein['modifications'] = trimmed_mods

        # Trim A3M MSAs
        if 'unpairedMsa' in protein and isinstance(protein['unpairedMsa'], str):
            protein['unpairedMsa'] = trim_a3m_string(protein['unpairedMsa'], protein_range)

        if 'pairedMsa' in protein and isinstance(protein['pairedMsa'], str):
            protein['pairedMsa'] = trim_a3m_string(protein['pairedMsa'], protein_range)

    with open(output_json_path, 'w') as f:
        json.dump(data, f, indent=2)

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
    elif type_glycan == "NAG(NAG(BMA(MAN)(MAN)))" or type_glycan == "M3":
        glycan_ccds = ['NAG', 'NAG', 'BMA', 'MAN', 'MAN']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']],
                       [[3,'O3'],[4, 'C1']],
                       [[3,'O6'],[5, 'C1']],
                       ]
    elif type_glycan == "NAG(NAG(BMA(MAN)(MAN(MAN)(MAN))))" or type_glycan == "M5":
        glycan_ccds = ['NAG', 'NAG', 'BMA', 'MAN', 'MAN', 'MAN', 'MAN']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']],
                       [[3,'O3'],[4, 'C1']],
                       [[3,'O6'],[5, 'C1']],
                       [[5, 'O3'], [6, 'C1']],
                       [[5, 'O6'], [7, 'C1']],
                       ]
    elif type_glycan == "NAG(NAG(BMA(MAN(MAN(MAN)))(MAN(MAN)(MAN(MAN)))))" or type_glycan == "M8":
        glycan_ccds = ['NAG', 'NAG', 'BMA', 'MAN', 'MAN', 'MAN', 'MAN', 'MAN', 'MAN', 'MAN']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']],
                       [[3,'O3'],[4, 'C1']],
                       [[3,'O6'],[5, 'C1']],
                       [[4, 'O2'], [6, 'C1']],
                       [[6, 'O2'], [7, 'C1']],
                       [[5, 'O3'], [8, 'C1']],
                       [[5, 'O6'], [9, 'C1']],
                       [[9, 'O2'], [10, 'C1']],
                       ]

    elif type_glycan == "NAG(FUC)(NAG(BMA(MAN)(MAN)))":
        glycan_ccds = ['NAG', 'NAG', 'BMA', 'MAN', 'MAN', 'FUC']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']],
                       [[3,'O3'],[4, 'C1']],
                       [[3,'O6'],[5, 'C1']],
                       [[1, 'O6'], [6, 'C1']],
                       ]
    elif type_glycan == "NAG(FUC)(NAG(BMA(MAN(NAG))(MAN(NAG))))" or type_glycan == "G0F":
        glycan_ccds = ['NAG', 'NAG', 'BMA', 'MAN', 'MAN', 'FUC', 'NAG', 'NAG']
        res_atom_id = 'ND2'
        glyc_connection_atom = [1, 'C1']
        inter_bonds = [[[1,'O4'],[2, 'C1']],
                       [[2,'O4'],[3, 'C1']],
                       [[3,'O3'],[4, 'C1']],
                       [[3,'O6'],[5, 'C1']],
                       [[1, 'O6'], [6, 'C1']],
                       [[4, 'O4'], [7, 'C1']],
                       [[5, 'O4'], [8, 'C1']],
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


def mask_template_region(json_path, chain_id, region_to_mask, template_num=0):
    """
    Deletes region that corresponds between the template sequence and query sequence, 1-based numbers
    :param json_path: Path to the AF3 JSON input file.
    :param chain_id: The chain ID of the protein which template should be modified.
    :param region_to_mask: tuple with two indices specifying the region to mask from the template sequence.
    :param template_num: Template number to musk (default is 0, 0-based numbers).
    """
    # TODO: rewrite function, it doesn't support multiple masking
    # Read the JSON file
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    # Find the protein chain with the specified chain_id
    chain_found = False
    for seq in af3_data.get("sequences", []):
        if "protein" in seq and seq["protein"].get("id") == chain_id:
            chain_found = True
            # Check if templates exist
            if "templates" not in seq["protein"] or len(seq["protein"]["templates"]) < template_num+1:
                raise ValueError(f"Template number {template_num} not found for chain {chain_id}")

            # Get the specific template
            template = seq["protein"]["templates"][template_num]
            query_indices = template["queryIndices"]
            template_indices = template["templateIndices"]

            # Validate region_to_mask
            start, end = region_to_mask
            if not (0 <= start <= end < len(template_indices)):
                raise ValueError(f"Invalid crop region {region_to_mask} for template of length {len(template_indices)}")

            # Create new lists excluding the cropped region
            new_query_indices = query_indices[:start] + query_indices[end:]
            new_template_indices = template_indices[:start] + template_indices[end:]

            # Update the template with new indices
            template["queryIndices"] = new_query_indices
            template["templateIndices"] = new_template_indices

            # Write the modified data back to the file
            with open(json_path, 'w') as f:
                json.dump(af3_data, f, indent=4)

            print(f"Cropped region {region_to_mask} from template {template_num} of chain {chain_id} in {json_path}")
            return

    if not chain_found:
        raise ValueError(f"Protein chain {chain_id} not found in {json_path}")

def fasta_to_json(fasta_path, output_json_path, model_seeds=(1, 2, 3, 4, 5), ids=("A", "B", "C", "D"),
                  dialect="alphafold3", version=1):
    with open(fasta_path, "r") as file:
        lines = file.readlines()

    if not lines or lines[0][0] != ">":
        raise ValueError("Invalid FASTA format: missing header line.")

    name = lines[0][1:].strip()
    sequence = "".join(line.strip() for line in lines[1:] if line.strip())

    if not sequence:
        raise ValueError("Invalid FASTA file: sequence is missing.")

    if sequence.count(">") > 0:
        raise ValueError("Multiple sequences detected in FASTA file. Only one is allowed.")

    json_data = {
        "name": name,
        "modelSeeds": model_seeds,
        "sequences": [
            {
                "protein": {
                    "sequence": sequence,
                    "id": ids
                }
            }
        ],
        "dialect": dialect,
        "version": version
    }

    with open(output_json_path, "w") as json_file:
        json.dump(json_data, json_file, indent=4)

    return json_data


# Example usage
if __name__ == "__main__":

    pure_json = "/g/kosinski/kgilep/flu_na_project/na_nc07/af3/input_json/na_nc07.json"

    structure_name = "t2cac4_optimized2.5"
    json_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/input_json/optimized2/{structure_name}.json"

    copy_input_json(pure_json, json_path, structure_name)

    na_chains = ["A", "B", "C", "D"]
    templates_path_dict = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized1/T2CAC4_full_chain_{id}.cif"
                           for id in na_chains}
    templates_path_dict_2 = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized1/T2CAC4_head_chain_{id}.cif"
                           for id in na_chains}
    templates_path_dict_3 = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized2/T2CAC4_ISOLDE_chain_A.cif"
                           for id in na_chains}
    paired_msa_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/msa/t2cac4_data_paired.a3m"
    unpaired_msa_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/msa/t2cac4_data_unpaired.a3m"
    query_range = (1, 468)
    region_to_mask_1 = (76,85)
    query_range_2 = (82,468)
    query_range_3 = (1, 468)
    region_to_mask_3 = (0, 80)

    # excluded 386 (not visible on the density map, can't see density under the Ab as well)
    glycosylation_dict = {42: 'G0F',
                          50: 'G0F',
                          58: 'G0F',
                          63: 'G0F',
                          68: 'G0F',
                          88: 'M3',
                          235: 'M3',
                          146: 'M3'}

    split_by_chains(json_path)
    change_input_json_version(json_path, 2)
    change_seeds(json_path, 10)

    for chain_id in na_chains:
        add_protein_template(json_path, chain_id, templates_path_dict[chain_id], query_range)
        mask_template_region(json_path, chain_id, region_to_mask_1, template_num=0)
        add_protein_template(json_path, chain_id, templates_path_dict_2[chain_id], query_range_2)
        add_protein_template(json_path, chain_id, templates_path_dict_3[chain_id], query_range_3)
        mask_template_region(json_path, chain_id, region_to_mask_3, template_num=2)
        for glycan_num, glycan_type in glycosylation_dict.items():
            add_glycan(json_path, chain_id, glycan_num, glycan_type)
        add_path_to_msa(json_path, chain_id, paired_msa_path, unpaired_msa_path)
