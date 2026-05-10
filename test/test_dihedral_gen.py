#!/usr/bin/env python
# coding: utf-8

import subprocess
import sys

# Run dihedral gen script to generate test output file
subprocess.run([sys.executable, '../generate_dihedral_csts.py', 'cdk2_17_1h1q.pdb', '17.params', '1h1q.params', 'test_atom_pair_csts', 'cstfiles/test_output_dihed_csts'])

# Compare test output to expected output
with open('cstfiles/test_output_dihed_csts') as f:
    test_file = f.read()

with open('expected_dihedral_csts_output') as f:
    expected_file = f.read()

if test_file == expected_file:
    print("Test passed - expected dihedral cst file and test-outputted file are equivalent.")
else:
    print("Test failed!")
    print("\nDifferences found - check cstfiles/test_output_dihed_csts against expected_dihedral_csts_output")
    sys.exit(1)