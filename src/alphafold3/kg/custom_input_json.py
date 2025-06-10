import json
import string
from pathlib import Path

try:
    from  alphafold3.common.folding_input import Input
except ImportError:
    print('Warning. Alphafold3 is not installed properly')


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

def clear_templates_for_chain(input_json_path, output_json_path, chain_id):
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

def trim_a3m_lines(a3m_lines, position_range):
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

    trimmed_lines = trim_a3m_lines(lines, position_range)

    with open(output_path, 'w') as f:
        f.writelines(trimmed_lines)

def trim_protein_chain(protein_range, chain_id, input_json_path, output_json_path):
    start, end = protein_range

    def trim_a3m_string(a3m_string, position_range):
        """Trim A3M block from a string using the same logic as A3M files."""
        lines = a3m_string.strip().splitlines(keepends=False)
        trimmed_lines = trim_a3m_lines([line + '\n' for line in lines], position_range)
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

def get_nglycan_positions(json_path, chain_id):
    """Finds the glycans' positions (1-based) that correspond to N
    in the Asn-X-Ser/Thr (X not Pro) consensus.

    :param json_path: Path to the AF3 JSON input file.
    :param chain_id: The chain ID to analyze for glycosylation sites.
    :return: Tuple of 1-based positions of potential N-glycosylation sites.
    """
    # Read the JSON file
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    # Find the sequence for the specified chain
    sequence = None
    for entry in af3_data["sequences"]:
        for key, value in entry.items():
            if isinstance(value, dict) and "id" in value:
                if isinstance(value["id"], list):
                    if chain_id in value["id"]:
                        sequence = value.get("sequence", "")
                        break
                elif value["id"] == chain_id:
                    sequence = value.get("sequence", "")
                    break
        if sequence:
            break

    if not sequence:
        return ()

    # List to store positions of N-glycosylation sites (1-based)
    positions = []

    # Check for N-X-S/T pattern (X is not Proline)
    for i in range(len(sequence) - 2):  # Need at least 3 residues for the motif
        if (sequence[i] == 'N' and
                sequence[i + 1] != 'P' and
                (sequence[i + 2] == 'S' or sequence[i + 2] == 'T')):
            positions.append(i + 1)  # 1-based indexing

    return tuple(positions)
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

def generate_indices(query_range: tuple[int]) -> tuple[list[int], list[int]]:
    """
    Generate query and template indices based on the given range.
    :param query_range: Tuple (start, end) for query indices.
    :return: queryIndices and templateIndices lists.
    """
    query_indices = list(range(query_range[0]-1, query_range[1]))
    template_indices = list(range(len(query_indices)))
    return query_indices, template_indices

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
    Deletes region that corresponds between the template sequence and query sequence, 1-based inclusive numbers.
    If the query index falls in the region_to_mask, remove both the query and corresponding template index.

    :param json_path: Path to the AF3 JSON input file.
    :param chain_id: The chain ID of the protein which template should be modified.
    :param region_to_mask: tuple with two indices (1-based, inclusive) specifying the region to mask.
    :param template_num: Template number to mask (default is 0, 0-based index).
    """
    # Read the JSON file
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    # Find the protein chain with the specified chain_id
    for seq in af3_data.get("sequences", []):
        if "protein" in seq and seq["protein"].get("id") == chain_id:
            templates = seq["protein"].get("templates", [])
            if len(templates) <= template_num:
                raise ValueError(f"Template number {template_num} not found for chain {chain_id}")

            template = templates[template_num]
            query_indices = template["queryIndices"]
            template_indices = template["templateIndices"]

            if len(query_indices) != len(template_indices):
                raise ValueError("Mismatch in length of queryIndices and templateIndices.")

            # 1-based inclusive masking range
            mask_start, mask_end = region_to_mask
            if mask_start > mask_end:
                raise ValueError("Invalid region_to_mask: start must be <= end")

            # Convert range to set of values to mask (0-based)
            mask_set = set(range(mask_start-1, mask_end))

            # Filter pairs
            new_query_indices = []
            new_template_indices = []

            for q_idx, t_idx in zip(query_indices, template_indices):
                if q_idx not in mask_set:
                    new_query_indices.append(q_idx)
                    new_template_indices.append(t_idx)

            # Update the template
            template["queryIndices"] = new_query_indices
            template["templateIndices"] = new_template_indices

            # Write back to JSON
            with open(json_path, 'w') as f:
                json.dump(af3_data, f, indent=4)

            print(f"Masked query positions {mask_start}-{mask_end} from template {template_num} of chain {chain_id} in {json_path}")
            return

    raise ValueError(f"Protein chain {chain_id} not found in {json_path}")


def fasta_to_homooligomer_json(
        fasta_path,
        output_json_path,
        structure_name=None,
        model_seeds=(1, 2, 3, 4, 5),
        ids=("A", "B", "C", "D"),
        dialect="alphafold3",
        version=2
):
    with open(fasta_path, "r") as file:
        lines = file.readlines()

    if not lines or lines[0][0] != ">":
        raise ValueError("Invalid FASTA format: missing header line.")

    if structure_name is not None:
        name = structure_name
    else:
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

    name_map = {
        "A0A0M4AQ34": "NB",
        "U5XIF8": "N2",
        "X2FPH5": "N3",
        "E4UH66": "N3",
        "D1LRJ8": "N4",
        "Q20VY9": "N5",
        "G0KN55": "N6",
        "A0A1U9GV95": "N7",
        "L8B2Z3": "N8",
        "E8Z0S3": "N9",
        "H6QM95": "N10",
        "U5N4D7": "N11",
    }

    fasta_dir = Path("/g/kosinski/kgilep/flu_na_project/na_variable/sequences")
    input_json_dir = Path("/g/kosinski/kgilep/flu_na_project/na_variable/af3/input_json")
    chain_ids = ["A", "B", "C", "D"]
    GLYCAN_TYPE = "M3"

    for fasta_path in fasta_dir.glob("*.fasta"):

        # Prepare raw JSON input from fasta
        prot_name = fasta_path.stem
        if prot_name in name_map:
            prot_name = f'{name_map[prot_name]}_{prot_name}'
        json_path = input_json_dir / (prot_name+".json")
        fasta_to_homooligomer_json(
                fasta_path,
                json_path,
                structure_name=prot_name,
                ids=chain_ids
        )
        # Add glycans
        glycan_positions = get_nglycan_positions(json_path, chain_ids[0])
        for glycan_pos in glycan_positions:
            for chain_id in chain_ids:
                add_glycan(json_path, chain_id, glycan_pos, GLYCAN_TYPE)
