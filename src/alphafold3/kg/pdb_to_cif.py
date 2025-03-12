import os
import sys
from Bio.PDB import PDBParser, MMCIFIO


def convert_pdb_to_cif(pdb_file):
    """Converts a PDB file to mmCIF format and saves it in the same directory."""

    if not os.path.isfile(pdb_file):
        print(f"Error: File '{pdb_file}' not found.")
        sys.exit(1)

    # Parse the structure
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("protein", pdb_file)

    # Generate output filename in the same directory
    output_cif = os.path.splitext(pdb_file)[0] + ".cif"

    # Save as mmCIF
    mmcif_io = MMCIFIO()
    mmcif_io.set_structure(structure)
    mmcif_io.save(output_cif)

    print(f"Conversion complete: {output_cif}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python pdb_to_cif.py <input.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    convert_pdb_to_cif(pdb_file)
