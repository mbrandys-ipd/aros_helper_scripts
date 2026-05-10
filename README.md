# Dihedral Constraint Generator

Generates dihedral constraints for use in an Alchemical Rosetta protocol. (Alchemical methods are used to predict the relative binding affinities of pairs of ligands in a lead optimization drug campaign).

Takes a protein-ligand complex PDB file where the two transformational ligands (ligand A and ligand B) are the last two residues, and outputs a dihedral constraint file for use in Rosetta. 

> **Note:** This tool is part of a larger alchemical Rosetta protocol currently in preparation for publication.

## Usage

```bash
python dihedral_gen.py input.pdb output.cst
```

## Testing

The `test/` folder contains example input files and an expected output. To verify your installation:

```bash
python test_dihedral_gen.py
```

This will generate a test output from the provided inputs, compare it to the expected output, and report whether the results match.

## Dependencies

- [PyRosetta](https://www.pyrosetta.org/)
- `contextlib` (standard library)
- `sys` (standard library)
