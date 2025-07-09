import json
import string
from pathlib import Path
import os
import shutil
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

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

    If paired_msa_path or unpaired_msa_path is empty or None, adds
    'pairedMsa' or 'unpairedMsa' fields with empty string instead of
    'pairedMsaPath' or 'unpairedMsaPath'.

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
                if paired_msa_path:
                    value["pairedMsaPath"] = paired_msa_path
                    value.pop("pairedMsa", None)  # Remove pairedMsa if present
                else:
                    value["pairedMsa"] = ""
                    value.pop("pairedMsaPath", None)  # Remove pairedMsaPath if present

                if unpaired_msa_path:
                    value["unpairedMsaPath"] = unpaired_msa_path
                    value.pop("unpairedMsa", None)  # Remove unpairedMsa if present
                else:
                    value["unpairedMsa"] = ""
                    value.pop("unpairedMsaPath", None)  # Remove unpairedMsaPath if present

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


def mutate_position(json_path: str, chain_id: str, mutation: str) -> None:
    """
    Introduces a mutation in an AlphaFold3 input JSON file for a specific chain.
    Only the first sequence in pairedMsa and unpairedMsa (matching the query sequence) is mutated.
    Handles MSA paths by reading the files and converting them to pairedMsa/unpairedMsa fields.

    Args:
        json_path: Path to the AF3 JSON input file.
        chain_id: Chain ID to mutate (e.g., 'A').
        mutation: Mutation in format 'positionAminoAcid' (e.g., '76N'), where position is 1-based.

    Raises:
        ValueError: If mutation format is invalid, chain not found, position out of range,
                   invalid amino acid, or first MSA sequence doesn't match query.
        FileNotFoundError: If MSA path files do not exist.
    """
    # Validate mutation format (e.g., '76N')
    mutation_pattern = r'^(\d+)([A-Z])$'
    match = re.match(mutation_pattern, mutation)
    if not match:
        raise ValueError(f"Invalid mutation format: '{mutation}'. Expected format: 'positionAminoAcid' (e.g., '76N')")

    position, new_aa = match.groups()
    position = int(position)

    # Validate amino acid
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    if new_aa not in valid_aa:
        raise ValueError(f"Invalid amino acid: '{new_aa}'. Must be one of {''.join(sorted(valid_aa))}")

    # Read JSON file
    json_path = Path(json_path)
    with open(json_path, 'r') as f:
        af3_data = json.load(f)

    # Find the protein chain
    protein_entry = None
    for entry in af3_data.get('sequences', []):
        if 'protein' in entry and entry['protein'].get('id') == chain_id:
            protein_entry = entry['protein']
            break
    if not protein_entry:
        raise ValueError(f"Chain ID '{chain_id}' not found in {json_path}")

    # Validate position and sequence
    sequence = protein_entry.get('sequence', '')
    if not sequence:
        raise ValueError(f"No sequence found for chain '{chain_id}'")
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} is out of range for sequence of length {len(sequence)}")

    # Apply mutation to sequence (0-based indexing internally)
    sequence_list = list(sequence)
    sequence_list[position - 1] = new_aa
    protein_entry['sequence'] = ''.join(sequence_list)

    # Process MSAs (from strings or paths)
    for msa_key, path_key in [('pairedMsa', 'pairedMsaPath'), ('unpairedMsa', 'unpairedMsaPath')]:
        msa = protein_entry.get(msa_key, '')
        msa_path = protein_entry.get(path_key)

        # If MSA path is present, read the file
        if msa_path:
            msa_path = Path(msa_path)
            if not msa_path.exists():
                raise FileNotFoundError(f"MSA file not found: {msa_path}")
            with open(msa_path, 'r') as f:
                msa = f.read()

        # Skip if no MSA content
        if not msa:
            continue

        # Split MSA into lines
        msa_lines = msa.strip().splitlines()
        if not msa_lines:
            continue

        # Process MSA sequences
        new_msa_lines = []
        is_first_sequence = True
        for line in msa_lines:
            if line.startswith('>'):
                new_msa_lines.append(line)
                is_first_sequence = False
                continue

            seq = line.strip()
            if is_first_sequence:
                # First sequence matches query, so use position directly
                if len(seq) != len(sequence):
                    raise ValueError(f"First sequence in {msa_key} does not match query sequence length")
                if seq != sequence:
                    raise ValueError(f"First sequence in {msa_key} does not match query sequence")
                seq_list = list(seq)
                seq_list[position - 1] = new_aa
                new_msa_lines.append(''.join(seq_list))
            else:
                # Keep subsequent sequences unchanged
                new_msa_lines.append(seq)

        # Update MSA in JSON and remove path if present
        protein_entry[msa_key] = '\n'.join(new_msa_lines) + '\n'
        if msa_path:
            protein_entry.pop(path_key, None)

    # Write updated JSON
    with open(json_path, 'w') as f:
        json.dump(af3_data, f, indent=4)

def download_uniprot_fasta(accession, output_dir):
    """
    Download a FASTA file for a single UniProt accession and save it to the output directory.

    Parameters:
        accession (str): UniProt accession string (e.g., "Q5BUA8")
        output_dir (str or Path): Path to the output directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        with urlopen(url) as response:
            fasta_data = response.read().decode("utf-8")
            fasta_path = output_dir / f"{accession}.fasta"
            fasta_path.write_text(fasta_data)
            print(f"Downloaded {accession} -> {fasta_path}")
    except HTTPError as e:
        print(f"HTTP error for {accession}: {e.code}")
    except URLError as e:
        print(f"URL error for {accession}: {e.reason}")

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

def relocate_and_symlink_jsons(base_dir, name_map):
    base_dir = Path(base_dir)
    data_dir = base_dir / "data"
    link_dir = base_dir / "link_data"

    data_dir.mkdir(exist_ok=True)
    link_dir.mkdir(exist_ok=True)

    for acc, shortname in name_map.items():
        pattern = f"{shortname.lower()}_{acc.lower()}"
        target_dir = base_dir / pattern

        if not target_dir.exists():
            print(f"Warning: {target_dir} does not exist.")
            continue

        json_files = list(target_dir.glob("*_data.json"))
        if not json_files:
            print(f"No *_data.json found in {target_dir}")
            continue

        for json_file in json_files:
            new_path = data_dir / json_file.name
            shutil.move(str(json_file), new_path)

            relative_target = os.path.relpath(new_path, start=link_dir)
            symlink_path = link_dir / json_file.name
            if symlink_path.exists() or symlink_path.is_symlink():
                symlink_path.unlink()
            symlink_path.symlink_to(relative_target)

            print(f"Moved: {json_file} -> {new_path}")
            print(f"Linked: {symlink_path} -> {relative_target}")

        try:
            target_dir.rmdir()
        except OSError:
            print(f"Could not remove {target_dir} (not empty or in use)")

# Example usage
if __name__ == "__main__":

    pure_json = "/g/kosinski/kgilep/flu_na_project/na_nc07/af3/input_json/na_nc07.json"

    structure_name = "t2cac4_optimized3.0"
    json_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/input_json/optimized3/{structure_name}.json"

    copy_input_json(pure_json, json_path, structure_name)

    na_chains = ["A", "B", "C", "D"]
    templates_path_dict = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized1/T2CAC4_full_chain_{id}.cif"
                           for id in na_chains}
    templates_path_dict_2 = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized1/T2CAC4_head_chain_{id}.cif"
                           for id in na_chains}
    templates_path_dict_3 = {id : f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/templates/optimized2/T2CAC4_ISOLDE_chain_A.cif"
                           for id in na_chains}
    paired_msa_path = f"/g/kosinski/kgilep/flu_na_project/na_nc07/af3/msa/T2CAC4sorted_na_group1_mafft_filtered_restored.a3m"
    unpaired_msa_path = ""
    query_range = (1, 468)
    region_to_mask_1 = (76,85)
    query_range_2 = (82,468)
    query_range_3 = (1, 468)
    region_to_mask_3 = (0, 76)

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

    for chain_id in na_chains:
        add_protein_template(json_path, chain_id, templates_path_dict[chain_id], query_range)
        mask_template_region(json_path, chain_id, region_to_mask_1, template_num=0)
        add_protein_template(json_path, chain_id, templates_path_dict_2[chain_id], query_range_2)
        add_protein_template(json_path, chain_id, templates_path_dict_3[chain_id], query_range_3)
        mask_template_region(json_path, chain_id, region_to_mask_3, template_num=2)
        for glycan_num, glycan_type in glycosylation_dict.items():
            add_glycan(json_path, chain_id, glycan_num, glycan_type)
        add_path_to_msa(json_path, chain_id, paired_msa_path, unpaired_msa_path)


## For prediction other NAs
# if __name__ == "__main__":
#
#     # name_map = {
#     #     "A0A0M4AQ34": "NB",
#     #     "U5XIF8": "N2",
#     #     "X2FPH5": "N3",
#     #     "E4UH66": "N3",
#     #     "D1LRJ8": "N4",
#     #     "Q20VY9": "N5",
#     #     "G0KN55": "N6",
#     #     "A0A1U9GV95": "N7",
#     #     "L8B2Z3": "N8",
#     #     "E8Z0S3": "N9",
#     #     "H6QM95": "N10",
#     #     "U5N4D7": "N11",
#     # }
#     # name_map = {
#     #     "A0A8K1EM63" : "N2"
#     # }
#     # name_map = {
#     #     "A3KE38" : "N3",
#     #     "C4LLY1" : "N4",
#     #     "N12" : "",
#     # }
#
#     name_map = {
#         "Q5BUA8": "N8_1",  # two PP, first subgroup of N8 NAs
#         "A0A023LRK0": "N8_main",  # most abundant, predicted was from this group
#         "A0A0B4V008": "N8_3",  # third subgroup
#         "E3JMK4": "N8_rare",  # only few such
#
#         "P03478": "N5_del",  # weird variant with the deletion underhead
#         "H8P788": "N5_other",  # different from the used
#         "G0KJY5": "N5_main",  # most abundant
#     }
#
#     name_map = {
#         # N6
#         "G7WV18": "N6_1_IKED",  # group 1 with IKED motif
#         "Q6XV50": "N6_1_del",  # group 1 with deletion
#         "X2G052": "N6_2_NKNE",  # group 2 with NKNE motif
#         "D6RUH4": "N6_2_del1",  # group 2 with deletion1
#         "F1BDM1": "N6_2_del2",  # group 2 with deletion2
#
#         # N9
#         "A0A0B4Q774": "N9_del",  # with deletion
#     }
#
#     name_map = {
#         # N3
#         "Q7TF26": "N3_1_del",  # group 1 with deletion
#
#         # N7
#         "A0A1S6GT84": "N7_1_main",  # group 1 (main)
#         "A0A0A7CIT3": "N7_1.5_del",  # deletion between group 1 and 2
#         "A0A024CQQ5": "N7_2",  # group 2
#         "P88837": "N7_3",  # small group 3
#     }
#
#     name_map = {
#         # N2
#         "V5SLV7": "N2_logo"
#     }
#
#     # fasta_dir = Path("/g/kosinski/kgilep/flu_na_project/na_variable/sequences")
#     # input_json_dir = Path("/g/kosinski/kgilep/flu_na_project/na_variable/af3/input_json")
#     # chain_ids = ["A", "B", "C", "D"]
#     # GLYCAN_TYPE = "M3"
#     # NUM_SEEDS = 10
#     #
#     # for uniprot, na_type in name_map.items():
#     #     download_uniprot_fasta(uniprot, fasta_dir)
#     #     prot_name = f'{na_type}_{uniprot}'
#     #     fasta_path = fasta_dir / f"{uniprot}.fasta"
#     #     json_path = input_json_dir / (prot_name+".json")
#     #     fasta_to_homooligomer_json(
#     #             fasta_path,
#     #             json_path,
#     #             structure_name=prot_name,
#     #             ids=chain_ids,
#     #             model_seeds=(tuple(range(NUM_SEEDS)))
#     #     )
#     #     # Add glycans
#     #     glycan_positions = get_nglycan_positions(json_path, chain_ids[0])
#     #     for glycan_pos in glycan_positions:
#     #         for chain_id in chain_ids:
#     #             add_glycan(json_path, chain_id, glycan_pos, GLYCAN_TYPE)
#
#     result_dir="/scratch/kgilep/flu_na_project/na_variable/af3"
#     relocate_and_symlink_jsons(result_dir, name_map)