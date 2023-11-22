#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
    """
    Determines the dihedral angles (phi,psi) for each residue in a protein and absolute precision error
    """

__author__ = "Hector Diaz"
__date__ = "08092021"

import os
import sys

import xlsxwriter

# from pdbtools2 import pdb_download
# from pdbtools2.helper import cmdline, geometry
import geometry


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
                if not resid_contents.has_key(line[13:16]):
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
                aux_dihedral.append(geometry.calcChiDihedrals(struct[residue][chi1atoms[res][0]], struct[residue][chi1atoms[res][1]],
                                                           struct[residue][chi1atoms[res][2]], struct[residue][chi1atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi2 and len(struct[residue]) > 4:
            try:
                aux_dihedral.append(geometry.calcChiDihedrals(struct[residue][chi2atoms[res][0]], struct[residue][chi2atoms[res][1]],
                                                           struct[residue][chi2atoms[res][2]], struct[residue][chi2atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi3 and len(struct[residue]) > 5:
            try:
                aux_dihedral.append(geometry.calcChiDihedrals(struct[residue][chi3atoms[res][0]], struct[residue][chi3atoms[res][1]],
                                                           struct[residue][chi3atoms[res][2]], struct[residue][chi3atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi4 and len(struct[residue]) > 6:
            try:
                aux_dihedral.append(geometry.calcChiDihedrals(struct[residue][chi4atoms[res][0]], struct[residue][chi4atoms[res][1]],
                                                           struct[residue][chi4atoms[res][2]], struct[residue][chi4atoms[res][3]]))
            except ValueError:
                err = "Dihedral calculation failed for %s\n" % residue_list[i]
                sys.stderr.write(err)
        else:
            aux_dihedral.append(0)

        if residue[0:3] in chi5 and len(struct[residue]) > 7:
            try:
                aux_dihedral.append(geometry.calcChiDihedrals(struct[residue][chi5atoms[res][0]], struct[residue][chi5atoms[res][1]],
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
                if not resid_contents.has_key(line[13:16]):
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
            dihedrals.append(geometry.calcDihedrals(CO[i - 1], N[i], CA[i], CO[i],
                                                    N[i + 1]))
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
    f = open(pdb_file, 'r')
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

def buscar_linea(rotamero_nativo, residuo, cadena):
    """
    Busca un rotamero en la lista de rotameros nativos
    """
    res = None
    for linea in rotamero_nativo:
        if linea[47:51]:
            if linea[46:47] == cadena:
                if linea[47:51] == residuo:
                    res = linea
    return res

def main():
    """
    Call if this is called from the command line.
    Calcula el porcentaje de residuos correctos segun la medida de precision absoluta basada en un angulo de 20
    ejecucion: python2 precision_absoluta.py list_pdb_origen list_pdb_faspr
    """

    # leer archivos con lista de ubicaciones de pdbs el primero debe contener a los originales
    # se espera que las listas esten ordenadas y tengan la misma cantidad de pdbs, en el caso que un pdb sea procesable
    # se debe colocar *** en la lista
    # todos los pdbs deben estar previamente estructurados de manera a tener la misma secuencia de residuos
    file_list = sys.argv[1:]
    pdbs_file_list = {}
    pdbs_file_key = []
    for pdb_list_file in file_list:
        pdbs_file_key.append(pdb_list_file)
        f = open(pdb_list_file, 'r')
        pdbs_file_list[pdb_list_file] = f.read().splitlines()
        f.close()

    print('#######################################################')
    print("///// Resultados de medida por Exactitud Absoluta /////")
    correct_angle = 20  # angulo tomado como correcto
    print_pdb = False  # imprimir pdbs en consola
    res_metodos = {} # diccionario para almacenar los resultados de cada metodo
    # creo un excel para almacenar los resultados
    workbook = xlsxwriter.Workbook('Resultados_presicion_absoluta_casp14_rota4_20.xlsx')
    worksheet = workbook.add_worksheet()
    row = 0
    col = 0
    worksheet.write(row, col, "Resulados de medida por Exactitud Absoluta")
    row += 1
    worksheet.write(row, col, "Angulo de tolerancia: %s" % correct_angle)
    row += 1

    # recorrer por listas para comparar la primera con cada una de las restantes
    # k es el identificador del archivo lista a comparar
    for k in range(1,len(pdbs_file_key)):
        worksheet.write(row, col, "Metodo")
        worksheet.write(row, col + 1, pdbs_file_key[k])
        row += 1
        worksheet.write(row, col, "PDB")
        worksheet.write(row, col + 1, "Residuos")
        worksheet.write(row, col + 2, "Evaluado X1")
        worksheet.write(row, col + 3, "Correcto X1")
        worksheet.write(row, col + 4, "Porcentaje X1")
        worksheet.write(row, col + 5, "Evaluado X1+2")
        worksheet.write(row, col + 6, "Correcto X1+2")
        worksheet.write(row, col + 7, "Porcentaje X1+2")
        row += 1

        # recorrer dentro de las listas para obtener los pdbs a comparar
        # h es el identificador del pdb dentro del archivo lista a comparar
        for h in range(len(pdbs_file_list[pdbs_file_key[k]])):
            # por cada lista se tiene un pdb
            pdb_file = pdbs_file_list[pdbs_file_key[0]][h]
            if pdb_file[:3] == "***":
                continue
            # almacenamos los rotameros del pdb original
            rotameros_nativo=torsionCHI(pdb_file)
            pdb_file = pdbs_file_list[pdbs_file_key[k]][h]
            # se agrega un filtro para saltar la comparacion si se encuentra una linea de *** en la lista de pdb
            if pdb_file[:3] == "***":
                continue
            # almacenamos los rotameros del pdb a comparar
            rotameros_metodo = torsionCHI(pdb_file)
            
            correctox1 = 0 # cantidad de correctos X1
            evaluadox1 = 0 # cantidad de residuos validos para evaluar x1
            correctox1x2 = 0 # cantidad de correctos X1X2
            evaluadox1x2 = 0 # cantidad de residuos validos para evaluar x1x2


            # agregar columnas
            rotameros_metodo[0] = rotameros_metodo[0][:-1]+"%10s%10s%15s%15s%10s%10s"   %("Ch1orig","Ch2orig","EvaluableX1","EvaluableX1+2","X1","X1+2\n")

            # out.append("%30s%4s \"%s\"%10.2F%10.2F%10.2F%10.2F%10.2F%10.2F%10.2F\n" % (short_pdb, labels[i][:3], labels[i][4:], phi_psi[i][0], phi_psi[i][1],dihedrals[i][0], dihedrals[i][1], dihedrals[i][2], dihedrals[i][3],dihedrals[i][4]))

            for j in range(1,len(rotameros_metodo)):
                # almacenamos los chi1 y chi2 originales y predichos
                # siendo el pdb native (rotameros_nativo) el pdb de origen
                # buscamos la linea que pertenece al mismo residuo en el pdb native para obtener los valores de chi1 y chi2
                linea_nativo = buscar_linea(rotameros_nativo,rotameros_metodo[j][47:51], rotameros_metodo[j][46:47])
                if not linea_nativo:
                    print ("Error: no se encontro el rotamero en el pdb nativo"+ rotameros_metodo[j][10:51])
                    continue
                ori_amino = linea_nativo[41:44]
                ori_chi1 = float(linea_nativo[72:82])
                ori_chi2 = float(linea_nativo[82:92])
                pred_amino = rotameros_metodo[j][41:44]
                pred_chi1= float(rotameros_metodo[j][72:82])
                pred_chi2 = float(rotameros_metodo[j][82:92])
                rotameros_metodo[j] = rotameros_metodo[j][:-1]+"%10.2F%10.2F"   %(ori_chi1,ori_chi2)

                #Calculo de diferencia entre angulos teniendo en cuenta el cuadrante en el que se encuentra.
                Sx1 = 360
                Sx1x2 = 360
                # consideraciones para evaluar x1 y x1x2
                if ori_chi1 != 0 and pred_chi1 != 0:
                    # Se resta el valor menor al mayor para obtener el valor absoluto de la diferencia
                    Ex1= max(ori_chi1,pred_chi1) - min(ori_chi1,pred_chi1)
                    # Se evalua si el valor absoluto de la diferencia es mayor a 180
                    Sx1= min(Ex1,360-Ex1)
                    evaluadox1 = evaluadox1 + 1 # aumentar la cantidad de residuos evaluados
                    rotameros_metodo[j] = rotameros_metodo[j] + "%15s" % ("1") # agregar el valor para la tabla
                    if ori_chi2 != 0 and pred_chi2 != 0:
                        # Consideraciones especiales segun ERAN EYAL (2004) por simetria de los aminoacidos para chi2
                        if pred_amino == 'ASN' or pred_amino == 'HIS' or pred_amino == 'ASP' or pred_amino == 'PHE' or pred_amino == 'TYR':
                            Ex1x2 = max(ori_chi2, pred_chi2) - min(ori_chi2, pred_chi2)
                            dif1 = min(Ex1x2, 360 - Ex1x2)
                            # debemos comparar la diferencia entre el valor de la predicion y su simetrico (rotado 180 grados)
                            if pred_chi2 > 0:
                                pred_chi2 = pred_chi2 - 180
                            else:
                                pred_chi2 = pred_chi2 + 180
                            Ex1x2 = max(ori_chi2, pred_chi2) - min(ori_chi2, pred_chi2)
                            dif2 = min(Ex1x2, 360 - Ex1x2)
                            # tomamos el que de menor diferencia
                            Sx1x2 = min(dif1,dif2)
                        else:
                            Ex1x2 = max(ori_chi2, pred_chi2) - min(ori_chi2, pred_chi2)
                            Sx1x2 = min(Ex1x2, 360 - Ex1x2)
                        evaluadox1x2 = evaluadox1x2 + 1 # aumentar la cantidad de residuos evaluados
                        rotameros_metodo[j] = rotameros_metodo[j] + "%15s" % ("1") # agregar el valor para la tabla
                    else:
                        rotameros_metodo[j] = rotameros_metodo[j] + "%15s" % ("0")  # agregar el valor para la tabla
                else:
                    rotameros_metodo[j] = rotameros_metodo[j] + "%15s%15s" % ("0","0")  # agregar el valor para la tabla
                if Sx1 <= correct_angle:
                    rotameros_metodo[j] = rotameros_metodo[j] + "%10s"%("1")
                    correctox1 = correctox1+1 # aumentar la cantidad de correctos para el calculo de la probabilidad
                    if Sx1x2 <= correct_angle:
                        rotameros_metodo[j] = rotameros_metodo[j] + "%10s"%("1\n")
                        correctox1x2 = correctox1x2 + 1
                    else:
                        rotameros_metodo[j] = rotameros_metodo[j] + "%10s"%("0\n")
                else:
                    rotameros_metodo[j] = rotameros_metodo[j] + "%10s%10s"%("0","0\n")

            prob_x1 = correctox1*100/evaluadox1
            prob_x1x2 = correctox1x2*100/evaluadox1x2
            # impresion de resultados
            estado_arte= pdbs_file_key[k]
            print("Estado del arte: %s" %(estado_arte))
            if pdbs_file_list[pdbs_file_key[0]][h].split('/'):
               if pdbs_file_list[pdbs_file_key[0]][h].split('/')[-1].split(".")[0]:
                   pdb_name = pdbs_file_list[pdbs_file_key[0]][h].split('/')[-1].split(".")[0]
               else:
                   pdb_name = pdbs_file_list[pdbs_file_key[0]][h].split('/')[-1]
            else:
                pdb_name = pdbs_file_list[pdbs_file_key[0]][h]
            print("PDB: %s" %(pdb_name))
            print("%20s%20s" % ("X1(%)", "X1+2(%)"))
            print( "%20.2F%20.2F" %(prob_x1, prob_x1x2))
            print("%10s%10.0F%20.0F" %("Evaluado  ", evaluadox1,evaluadox1x2))
            print("%10s%10.0F%20.0F" %("Correcto  ",correctox1,correctox1x2))

            #Sumatorias por estado del arte
            if not estado_arte in res_metodos:
                res_metodos[estado_arte]={}
                res_metodos[estado_arte]={
                    'metodo': estado_arte,
                    'residuos': len(rotameros_metodo)-1,
                    'evaluadox1': evaluadox1,
                    'evaluadox1x2': evaluadox1x2,
                    'correctox1': correctox1,
                    'correctox1x2': correctox1x2,
                }
            else:
                res_metodos[estado_arte] = {
                    'metodo': estado_arte,
                    'residuos': res_metodos[estado_arte]['residuos']+len(rotameros_metodo) - 1,
                    'evaluadox1': res_metodos[estado_arte]['evaluadox1']+evaluadox1,
                    'evaluadox1x2': res_metodos[estado_arte]['evaluadox1x2']+evaluadox1x2,
                    'correctox1': res_metodos[estado_arte]['correctox1']+correctox1,
                    'correctox1x2': res_metodos[estado_arte]['correctox1x2']+correctox1x2,
                }
            # imprimir en consola
            if print_pdb:
                print()
                print("".join(rotameros_metodo))
            print_to_excel(rotameros_metodo, pdb_name, estado_arte)

            worksheet.write(row, col, pdb_name)
            worksheet.write(row, col + 1, len(rotameros_metodo) - 1)
            worksheet.write(row, col + 2, evaluadox1)
            worksheet.write(row, col + 3, correctox1)
            worksheet.write(row, col + 4, prob_x1)
            worksheet.write(row, col + 5, evaluadox1x2)
            worksheet.write(row, col + 6, correctox1x2)
            worksheet.write(row, col + 7, prob_x1x2)
            row += 1

    row += 1
    worksheet.write(row, col, "Totales")
    row += 1
    worksheet.write(row, col, "Metodo")
    worksheet.write(row, col + 1, "Residuos")
    worksheet.write(row, col + 2, "Evaluado X1")
    worksheet.write(row, col + 3, "Correcto X1")
    worksheet.write(row, col + 4, "Porcentaje X1")
    worksheet.write(row, col + 5, "Evaluado X1+2")
    worksheet.write(row, col + 6, "Correcto X1+2")
    worksheet.write(row, col + 7, "Porcentaje X1+2")
    row += 1
    # imprimir res_metodos
    for est_arte in res_metodos:
        print('Metodo: %s' %(res_metodos[est_arte]['metodo']))
        print('%15s%15s%15s%15s%25s' % ('Residuos', 'Evaluado X1', 'Evaluado X1+2', 'Correcto X1', 'Correcto X1+2'))
        print('%15.0F%15.0F%15.0F%15.0F%15.0F' % (res_metodos[est_arte]['residuos'], res_metodos[est_arte]['evaluadox1'], res_metodos[est_arte]['evaluadox1x2'], res_metodos[est_arte]['correctox1'], res_metodos[est_arte]['correctox1x2']))
        p1 = float(res_metodos[est_arte]['correctox1'])*100/float(res_metodos[est_arte]['evaluadox1'])
        p2 = float(res_metodos[est_arte]['correctox1x2'])*100/float(res_metodos[est_arte]['evaluadox1x2'])
        # print(p1)
        # print(p2)
        print('%15s%15.2F%17s%15.2F' % ('Probabilidad X1', p1, 'Probabilidad X1+2', p2))
        print('#############################################')
        print()
        worksheet.write(row, col, res_metodos[est_arte]['metodo'])
        worksheet.write(row, col + 1, res_metodos[est_arte]['residuos'])
        worksheet.write(row, col + 2, res_metodos[est_arte]['evaluadox1'])
        worksheet.write(row, col + 3, res_metodos[est_arte]['correctox1'])
        worksheet.write(row, col + 4, p1)
        worksheet.write(row, col + 5, res_metodos[est_arte]['evaluadox1x2'])
        worksheet.write(row, col + 6, res_metodos[est_arte]['correctox1x2'])
        worksheet.write(row, col + 7, p2)
        row += 1
    workbook.close()

#imprimir excel desde una lista de str formateando por celdas
def print_to_excel(rotameros, pdb_name, estado_arte):

    # crear directorio estado_arte si no existe
    if not os.path.exists('resultados/casp14/'+estado_arte):
        os.makedirs('resultados/casp14/'+estado_arte)
    xls_name= 'resultados/casp14/'+estado_arte+'/pa_'+pdb_name+"_"+estado_arte+'.xlsx'
    workbook = xlsxwriter.Workbook(xls_name)
    worksheet = workbook.add_worksheet()
    row = 0

    for line in (rotameros):

        col = 0
        #columna 1 numero
        if row == 0:
            celda = "N"
        else:
            celda = line[0:10]
        worksheet.write(row, col, celda)
        # columna 2 pdb
        celda = line[11:40]
        worksheet.write(row, col + 1, celda)
        # columna 3 aa
        celda = line[41:44]
        worksheet.write(row, col + 2, celda)
        # columna 4 res
        celda = line[45:52]
        worksheet.write(row, col + 3, celda)
        # columna 5 al 13 : phi, psi, chi1, chi2, chi3, chi4, chi5, chi6, chi1orig, chi2orig, X1, X1+2
        col=4
        for celda in line[53:143].split():
            worksheet.write(row, col, celda)
            col += 1
        # columna 14 al 17 EvaluadoX1, EvaluadoX1+2, X1, X1+2
        for celda in line[143:].split():
            if row>0:
                celda=int(celda)
            worksheet.write(row, col, celda)
            col += 1
        row += 1
    workbook.close()


if __name__ == "__main__":
    main()

