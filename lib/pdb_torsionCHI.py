#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
    """
    Determines the dihedral angles (phi,psi) for each residue in a protein and absolute precision error
    """

__author__ = "Carlos Bobadilla"
__date__ = "070601"

import os
import sys
from .geometry import calcChiDihedrals, calcDihedrals

def pdbTorsionSideChain(pdb):
    """
    Calculate the side-chain torsion angles for a pdb file.
    """

    residue_list = []
    #creating dictionary
    struct = {}



    # nor ALA neither GLY are included on any chi so values set to NA

    resid_contents = {}
    current_residue = None
    to_take = ["N  ", "CA ", "CB ", "CD ", "CD1", "CE ", "CG ", "CG1", "CZ ",
                "OG ", "OG1", "OD1", "OE1", "ND1", "NE ", "NH1", "NZ ", "SD ", "SG "]

    for line in pdb:
        if line[0:4] == "ATOM" or (line[0:6] == "HETATM" and line[17:20] == "MSE"):

            if line[13:16] in to_take:

                # First residue
                if current_residue == None:
                    current_residue = line[17:26]

                # If we're switching to a new residue, record the previously
                # recorded one.
                if current_residue != line[17:26]:

                    try:
                        aux = {}
                        for key in resid_contents:
                            positions = []
                            for i in range(3):
                                positions.append(float(resid_contents[key][30+8*i:39+8*i]))
                            aux[key] = positions

                        struct[current_residue] = aux
                        residue_list.append(current_residue)

                    except KeyError:
                        err = "Residue %s has missing atoms: skipping.\n" % current_residue
                        sys.stderr.write(err)

                    # Reset resid contents dictionary
                    current_residue = line[17:26]
                    resid_contents = {}

                # Now record N, C, and CA entries.  Take only a unique one from
                # each residue to deal with multiple conformations etc.
                if not line[13:16] in resid_contents:
                    resid_contents[line[13:16]] = line
                else:
                    err = "Warning: %s has repeated atoms!\n" % current_residue
                    sys.stderr.write(err)

    # Record the last residue
    try:
        aux = {}
        for key in resid_contents:
            positions = []
            for i in range(3):
                positions.append(float(resid_contents[key][30 + 8 * i:39 + 8 * i]))
            aux[key] = positions

        struct[current_residue] = aux
        residue_list.append(current_residue)

    except KeyError:
        err = "Residue %s has missing atoms: skipping.\n" % current_residue
        sys.stderr.write(err)

    # Calculate phi and psi for each residue.  If the calculation fails, write
    # that to standard error and move on.
    labels = []
    dihedrals = []
#chi1
    chi1atoms = {}
    chi1 = ["ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
            "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    chi1e = ["CYS", "ILE", "SER", "THR", "VAL"]

    for value in chi1:
        if value not in chi1e:
            atoms = ["N  ", "CA ", "CB ", "CG "]
            chi1atoms[value] = atoms
        elif value == 'CYS':
            atoms = ["N  ", "CA ", "CB ", "SG "]
            chi1atoms[value] = atoms
        elif value == 'ILE' or value == 'VAL':
            atoms = ["N  ", "CA ", "CB ", "CG1"]
            chi1atoms[value] = atoms
        elif value == 'SER':
            atoms = ["N  ", "CA ", "CB ", "OG "]
            chi1atoms[value] = atoms
        elif value == 'THR':
            atoms = ["N  ", "CA ", "CB ", "OG1"]
            chi1atoms[value] = atoms
#chi2
    chi2atoms = {}
    chi2 = ["ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "TRP", "TYR"]
    chi2e = ["ASN", "ASP", "HIS", "ILE", "LEU", "MET", "PHE", "TRP", "TYR"]

    for value in chi2:
        if value not in chi2e:
            atoms = ["CA ", "CB ", "CG ", "CD "]
            chi2atoms[value] = atoms
        elif value == 'ASN' or value == 'ASP':
            atoms = ["CA ", "CB ", "CG ", "OD1"]
            chi2atoms[value] = atoms
        elif value == 'HIS':
            atoms = ["CA ", "CB ", "CG ", "ND1"]
            chi2atoms[value] = atoms
        elif value == 'LEU' or value == 'PHE' or value == 'TRP' or value == 'TYR':
            atoms = ["CA ", "CB ", "CG ", "CD1"]
            chi2atoms[value] = atoms
        elif value == 'MET':
            atoms = ["CA ", "CB ", "CG ", "SD "]
            chi2atoms[value] = atoms
        elif value == 'ILE':
            atoms = ["CA ", "CB ", "CG1", "CD1"]
            chi2atoms[value] = atoms

#chi3

    chi3 = ["ARG", "GLN", "GLU", "LYS", "MET"]
    chi3atoms = {'ARG': ['CB ', 'CG ', 'CD ', 'NE '], 'GLN': ['CB ', 'CG ', 'CD ', 'OE1'],
                  'GLU': ['CB ', 'CG ', 'CD ', 'OE1'], 'LYS': ['CB ', 'CG ', 'CD ', 'CE '],
                  'MET':['CB ', 'CG ', 'SD ', 'CE ']}

#chi4

    chi4 = ["ARG", "LYS"]
    chi4atoms = {'ARG': ['CG ', 'CD ', 'NE ', 'CZ '], 'LYS': ['CG ', 'CD ', 'CE ', 'NZ ']}

#chi5
    chi5 = ["ARG"]
    chi5atoms = {'ARG': ['CD ', 'NE ', 'CZ ', 'NH1']}

#aux
    aux_dihedral = []
    i = 0
    cont = 0    #to count until len(residue_list) so we can break before the last atom.
    for residue in residue_list:    #
        cont = cont + 1
        if cont >= len(residue_list):
            break   #if we are on the last atom

        if cont == 1 and len(struct[residue]) < 4:
            continue   #if we are on the first atom and does not have a complete set for calcChiDihedrals

        res = residue[0:3]
        if residue[0:3] in chi1 and len(struct[residue]) >= 4:
            try:
                aux_dihedral.append(calcChiDihedrals(struct[residue][chi1atoms[res][0]], struct[residue][chi1atoms[res][1]],
                                              struct[residue][chi1atoms[res][2]], struct[residue][chi1atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi2 and len(struct[residue]) > 4:
            try:
                aux_dihedral.append(calcChiDihedrals(struct[residue][chi2atoms[res][0]], struct[residue][chi2atoms[res][1]],
                                              struct[residue][chi2atoms[res][2]], struct[residue][chi2atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi3 and len(struct[residue]) > 5:
            try:
                aux_dihedral.append(calcChiDihedrals(struct[residue][chi3atoms[res][0]], struct[residue][chi3atoms[res][1]],
                                              struct[residue][chi3atoms[res][2]], struct[residue][chi3atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi4 and len(struct[residue]) > 6:
            try:
                aux_dihedral.append(calcChiDihedrals(struct[residue][chi4atoms[res][0]], struct[residue][chi4atoms[res][1]],
                                              struct[residue][chi4atoms[res][2]], struct[residue][chi4atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi5 and len(struct[residue]) > 7:
            try:
                aux_dihedral.append(calcChiDihedrals(struct[residue][chi5atoms[res][0]], struct[residue][chi5atoms[res][1]],
                                              struct[residue][chi5atoms[res][2]], struct[residue][chi5atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        labels.append(residue_list[i])
        dihedrals.append(aux_dihedral)
        aux_dihedral = []
        i = i + 1


    return dihedrals, labels


def pdbTorsion(pdb):
    """
    Calculate the backbone torsion angles for a pdb file.
    """

    residue_list = []
    N = []
    CO = []
    CA = []

    resid_contents = {}
    current_residue = None
    to_take = ["N  ", "CA ", "C  "]
    for line in pdb:
        if line[0:4] == "ATOM" or (line[0:6] == "HETATM" and line[17:20] == "MSE"):

            if line[13:16] in to_take:

                # First residue
                if current_residue == None:
                    current_residue = line[17:26]

                # If we're switching to a new residue, record the previously
                # recorded one.
                if current_residue != line[17:26]:

                    try:
                        N.append([float(resid_contents["N  "][30 + 8 * i:39 + 8 * i])
                                  for i in range(3)])
                        CO.append([float(resid_contents["C  "][30 + 8 * i:39 + 8 * i])
                                   for i in range(3)])
                        CA.append([float(resid_contents["CA "][30 + 8 * i:39 + 8 * i])
                                   for i in range(3)])
                        residue_list.append(current_residue)

                    except KeyError:
                        err = "Residue %s has missing atoms: skipping.\n" % current_residue
                        sys.stderr.write(err)

                    # Reset resid contents dictionary
                    current_residue = line[17:26]
                    resid_contents = {}

                # Now record N, C, and CA entries.  Take only a unique one from
                # each residue to deal with multiple conformations etc.
                if not line[13:16] in resid_contents:
                    resid_contents[line[13:16]] = line
                else:
                    err = "Warning: %s has repeated atoms!\n" % current_residue
                    sys.stderr.write(err)

    # Record the last residue
    try:
        N.append([float(resid_contents["N  "][30 + 8 * i:39 + 8 * i])
                  for i in range(3)])
        CO.append([float(resid_contents["C  "][30 + 8 * i:39 + 8 * i])
                   for i in range(3)])
        CA.append([float(resid_contents["CA "][30 + 8 * i:39 + 8 * i])
                   for i in range(3)])
        residue_list.append(current_residue)

    except KeyError:
        err = "Residue %s has missing atoms: skipping.\n" % current_residue
        sys.stderr.write(err)

    # Calculate phi and psi for each residue.  If the calculation fails, write
    # that to standard error and move on.
    labels = []
    dihedrals = []
    for i in range(0, len(residue_list) - 1):
        try:
            dihedrals.append(calcDihedrals(CO[i - 1], N[i], CA[i], CO[i], N[i + 1]))
            labels.append(residue_list[i])
        except ValueError:
            err = "Dihedral calculation failed for %s\n" % residue_list[i]
            sys.stderr.write(err)

    dihedrals.append([0, 0])
    return dihedrals, labels

def torsionCHI(pdb_file):
    """
    Calcula los rotameros chi pertenecientes a cada residuo
    """
    out = []
    # Read in input file
    f = open(os.path.abspath('.')+"/"+pdb_file, 'r')
    pdb = f.readlines()
    f.close()

    # Calculate torsion angles and secondary structure
    dihedrals, labels = pdbTorsionSideChain(pdb)
    phi_psi, list_amino = pdbTorsion(pdb)

    # Print out results in pretty fashion
    short_pdb = os.path.split(pdb_file)[-1][:-4]
    for i in range(len(dihedrals)):
        out.append("%30s%4s \"%s\"%10.2F%10.2F%10.2F%10.2F%10.2F%10.2F%10.2F\n" % \
                   (short_pdb, labels[i][:3], labels[i][4:], phi_psi[i][0], phi_psi[i][1],
                    dihedrals[i][0], dihedrals[i][1], dihedrals[i][2], dihedrals[i][3],
                    dihedrals[i][4] if len(dihedrals[i]) == 5 else 0))


    out = ["%10i%s" % (i, x) for i, x in enumerate(out)]

    header = "%10s%30s%4s%8s%10s%10s%10s%10s%10s%10s%10s\n" % (
    " ", "pdb", "aa", "res", "phi", "psi", "chi1", "chi2", "chi3", "chi4", "chi5")
    out.insert(0, header)
    return out

