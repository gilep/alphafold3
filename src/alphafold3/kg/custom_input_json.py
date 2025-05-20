import json

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

            # Convert range to set of values to mask (1-based)
            mask_set = set(range(mask_start, mask_end + 1))

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

    ##################################################
    # BASIC EXAMPLE
    ##################################################
    json_raw = 'test/af3_input.json'
    json_path = 'test/af3_input_custom.json'  # file will be modified!!!
    template_path = '/user/test/template.cif' # should be monomeric, should be full path
    structure_name = "test_name"

    chain_ids = ["A"]
    query_range = (-40, 96)               # Range of the query residues if the second residue of the query corresponds
                                        # to the first residue of the template sequence.

    region_to_mask = (-40,40)             # Region in query that will not use the data from the template

    copy_input_json(json_raw, json_path, structure_name)
    split_by_chains(json_path)

    for chain_id in chain_ids:
        add_protein_template(json_path, chain_id, template_path, query_range)
        mask_template_region(json_path, chain_id, region_to_mask, template_num=0)

    ##################################################
    # FULL EXAMPLE WITH HOMOOLIGOMER AND MULTIPLE TEMPLATES
    ##################################################
    fasta_path = 'test/fasta.fasta'                      # input fasta wiht one sequnce!
    json_path = 'test/af3_input_homooligomer.json'       # file will be modified!!!
    template_path_1 = '/home/test/template.cif'          # should be monomeric
    template_path_2 = '/home/test/template2.cif'         # should be monomeric

    structure_name = "test"
    chain_ids = ["A", "B", "C", "D"]

    query_range_1 = (2, 96)
    region_to_mask_1 = (1,40)

    query_range_2 = (1, 96)
    region_to_mask_2 = (41,96)

    fasta_to_homooligomer_json(
        fasta_path,
        json_path,
        structure_name=structure_name,
        ids=chain_ids,
        model_seeds=(1,2,3)
    )

    ###############
    # run AF3 with --norun_inference with
    # to get test_data.json
    # replace af3_input_homooligomer.json with test_data.json
    # `cp path/to/test_data.json path/to/af3_input_homooligomer.json`
    ###############

    split_by_chains(json_path)

    for chain_id in chain_ids:
        # delete templates found by AF3
        clear_templates_for_chain(json_path, json_path, chain_id)

        add_protein_template(json_path, chain_id, template_path_1, query_range_1)
        mask_template_region(json_path, chain_id, region_to_mask_1, template_num=0)

        add_protein_template(json_path, chain_id, template_path_2, query_range_2)
        mask_template_region(json_path, chain_id, region_to_mask_2, template_num=1)

        # optional, it is recommended to let Alphafold3 generate MSA
        # add_path_to_msa(json_path, chain_id, paired_msa_path=paired_msa_path, unpaired_msa_path=unpaired_msa_path)
