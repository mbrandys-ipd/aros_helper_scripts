#!/usr/bin/env python
# coding: utf-8

from pyrosetta import *
from pyrosetta.rosetta import core
import os
import json
import numpy as np
import contextlib
import sys

'''
This script defines dihedral constraints for use in Alchemical Rosetta protocols; dihedrals are excluded where atoms already have an atom pair constraint defined. This script aims to find all possible dihedrals to constrain, but limiting the search pool based on the params-file specified dihedrals (CHIs). Atoms participating in a triple-bond are also excluded.
'''

# Hard-coded configurations
std_dev = "1.0" # standard deviation of constraint definition

def get_triple_bond_atoms(fn):
    '''
    Parse params file (fn) for any atoms participating in a triple bond.
    Returns a list of atom names so they can be excluded from dihedral constraints.
    '''
    with open(fn) as f:
        lines = f.readlines()

    exclude_these = []
    for line in lines:
        line = line.strip().split()
        if line[0] == "BOND_TYPE" and line[3] == "3":
            exclude_these.append(line[1])
            exclude_these.append(line[2])
            
    # clean up any redundant names found:
    uniques = list(set(exclude_these))
    
    return uniques

def get_chis(fn):
    '''
    Parse params file (fn) for 'CHI' dihedral definitions.
    For a dihedral defined as X1-atom1-atom2-X2, returns a list of [atom1, atom2], representing the core bond of each dihedral.
    '''
    with open(fn) as f:
        lines = f.readlines()

    chis = []
    for line in lines:
        line = line.strip().split()
        if line[0] == "CHI":
            chis.append(line[3:5])
    return chis

def get_bonds(fn):
    '''
    Parse params file (fn) for bond definitions using the 'BOND_TYPE' lines.
    Returns a list of lists (bonds) that has the atom names of bonded pairs.
    '''
    with open(fn) as f:
        lines = f.readlines()

    bonds = []
    for line in lines:
        line = line.strip().split()
        if line[0] == "BOND_TYPE":
            bonds.append(line[1:3])
    return bonds

def get_atomX_set(atom, excluded_atom, bonds):
    '''
    Given a dihedral defined like:
    X1 - atom1 - atom2 - X2
    this function returns a list (X_set) of all atoms that could cap the dihedral depending on what's passed in as 'atom' and 'excluded_atom':
    - If atom=atom1, excluded_atom=atom2, and we return candidates for X1 in X_set
    - If atom=atom2, excluded_atom=atom1, and we return candidates for X2 in X_set
    '''
    X_set = []
    for bond_pair in bonds:
        if atom in bond_pair:
            a = [x for x in bond_pair if ((x != atom) and (x != excluded_atom))]
            if len(a) > 0:
                X_set.append(a[0])
    return X_set

def get_dihedral(pose,resnum,atom_list):
    """Calculate dihedral angle for four atoms in radians."""
    
    # Get atom coordinates - uses pyrosetta commands
    p0 = np.array([pose.residue(resnum).atom(atom_list[0]).xyz()[0],pose.residue(resnum).atom(atom_list[0]).xyz()[1],pose.residue(resnum).atom(atom_list[0]).xyz()[2]]) 
    p1 = np.array([pose.residue(resnum).atom(atom_list[1]).xyz()[0],pose.residue(resnum).atom(atom_list[1]).xyz()[1],pose.residue(resnum).atom(atom_list[1]).xyz()[2]]) 
    p2 = np.array([pose.residue(resnum).atom(atom_list[2]).xyz()[0],pose.residue(resnum).atom(atom_list[2]).xyz()[1],pose.residue(resnum).atom(atom_list[2]).xyz()[2]]) 
    p3 = np.array([pose.residue(resnum).atom(atom_list[3]).xyz()[0],pose.residue(resnum).atom(atom_list[3]).xyz()[1],pose.residue(resnum).atom(atom_list[3]).xyz()[2]]) 

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.arctan2(y, x)

def check_for_apc(apc_list,chi):
    """
    Check if any atoms in the proposed dihedral already have atom pair constraints.
    Returns True if dihedral can be created (no overlap with APC atoms).
    """
    create = False

    chi_atoms = chi.split()
    if not any(apc in apc_list for apc in chi_atoms):
        create = True
    
    return create

def parse_apc_file(apc_fn, ligA_resnum, ligB_resnum):
    """
    Parse atom pair constraint file (apc_fn) to get lists of constrained atoms for each ligand.
    Returns lists of atom names with APC assigned.
    """
    with open(apc_fn) as f:
        lines = f.readlines()
    
    ligA_apc = []
    ligB_apc = []
    
    for line in lines:
        parts = line.strip().split()
        resnum = int(parts[4])
        atom = parts[3]
        
        if resnum == ligA_resnum and atom not in ligA_apc:
            ligA_apc.append(atom)
        if resnum == ligB_resnum and atom not in ligB_apc:
            ligB_apc.append(atom)
    
    return ligA_apc, ligB_apc

def gen_dihedral_csts(pose, params_file, lig_apc, lig_resnum, std_dev):
    '''
    Generate the dihedral constraints and save to a list that is returned.
    '''
    cst_list = []
    
    chis = get_chis(params_file)
    bonds = get_bonds(params_file)

    for chi in chis:
        atom1 = chi[0]
        atom2 = chi[1]
        X1s = get_atomX_set(atom1,atom2,bonds)
        X2s = get_atomX_set(atom2,atom1,bonds) 
        if X1s != None and X2s != None:
            for x1 in X1s:
                for x2 in X2s:
                    # checks if we can even make the dihedral based on if any atom is already APC-constrained
                    create_dihed = check_for_apc(lig_apc,f'{x1} {atom1} {atom2} {x2}')
                    if create_dihed:
                        # even if we can make the dihedral we need to exclude any that might contain our triple bonded atoms:
                        lig_tbond_atoms = get_triple_bond_atoms(params_file)
                        if any(item in [x1,atom1,atom2,x2] for item in lig_tbond_atoms):
                            continue
                        else:
                            # finally, calculate the dihedral and save the line
                            tor = get_dihedral(pose, lig_resnum, [x1,atom1,atom2,x2])
                            cst_list.append(f"Dihedral {x1} {lig_resnum} {atom1} {lig_resnum} {atom2} {lig_resnum} {x2} {lig_resnum} CIRCULARHARMONIC {tor:.2f} {std_dev}\n")
                    else:
                        continue
    return cst_list

def main():
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Usage: python generate_dihedral_csts.py <pdb_file> <ligandA_params_file> <ligandB_params_file> <apc_cstfile> [output_filename]")
        sys.exit(1)
    
    # path to pdb; expects a complex with ligandA and ligandB as the last two residues
    posePDB = sys.argv[1]
    # paths to rosetta params files for each ligand
    paramsFileA = sys.argv[2]
    paramsFileB = sys.argv[3]
    # path to atom pair constraint file generated ahead of time
    apc_fn = sys.argv[4]
    
    cstfile_fn = sys.argv[5] if len(sys.argv) == 6 else f"dihed_cst_stdev_{std_dev}"
    
    cstdir = os.path.join("cstfiles")
    os.makedirs(cstdir, exist_ok=True)
    
    # this will be the path to the outputted constraint file
    cstfile_path = os.path.join(cstdir, cstfile_fn)

    overwrite = False # this will overwrite the cstfile if it exists already
    if overwrite or not os.path.exists(cstfile_path):

        if not os.path.exists(apc_fn):
            raise FileNotFoundError(f"Atom pair cst file was not found. You entered: {apc_fn}.")

        # this block initializes pyrosetta and the starting pose
        args = f"-gen_potential -beta_cart -extra_res_fa {paramsFileA} {paramsFileB} -mute all"
        # this next section blocks anything printed to console from pyrosetta just to allow
        # this script to be integrated into a wrapper
        # where the wrapper captures console outputs
        with open(os.devnull, 'w') as fnull:
            with contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
                pyrosetta.init(args)
                
        # stores pose from pdbfile
        pose = pyrosetta.pose_from_pdb(posePDB)
        
        # sets up residue numbering of ligands
        ligA_res = pose.residue(pose.total_residue() - 1)
        ligA_resnum = ligA_res.seqpos()
        ligB_res = pose.residue(pose.total_residue())
        ligB_resnum = ligB_res.seqpos()
        
        # gets the atoms involved in atom pair constraints
        ligA_apc, ligB_apc = parse_apc_file(apc_fn, ligA_resnum, ligB_resnum)

        # collects the dihedral constraint lines for each ligand before combining
        lig_a_cst_set = gen_dihedral_csts(pose, paramsFileA, ligA_apc, ligA_resnum, std_dev)
        lig_b_cst_set = gen_dihedral_csts(pose, paramsFileB, ligB_apc, ligB_resnum, std_dev)
        cst_set = lig_a_cst_set + lig_b_cst_set

        # if no constraints could be assigned, tell the user
        if not cst_set:
            raise RuntimeError("No dihedral constraints met criteria - output would be empty. Exiting.")
        # otherwise write the file:
        with open(cstfile_path, 'w') as outfh:
            outfh.writelines(cst_set)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        error_output = {"error": str(e)}
        print(json.dumps(error_output))
        sys.exit(1)