#!/usr/bin/env python

# Copyright 2023, Hector Diaz
# This program is distributed under General Public License v. 3.
# COPYING for a copy of the license.

__description__ = \
    """
    Determine the RMSD (Root-mean-square deviation)
    """

__author__ = "Hector Diaz"
__date__ = "22022023"

import os
from Bio.PDB import PDBParser, Superimposer


def calcular_rmsd_ecuacion(atoms_structure1, atoms_structure2):
    """RMSD = raiz cuadrada de la suma de los cuadrados de las diferencias de las coordenadas
       dividido entre el numero de atomos
       RMSD = sqrt(1/N * sum((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2))
    """
    N = 0
    suma = 0
    for i in range(len(atoms_structure1)):
        for atom1, atom2 in zip(atoms_structure1[i], atoms_structure2[i]):
            x1, y1, z1 = atom1.get_coord()
            x2, y2, z2 = atom2.get_coord()
            suma += ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            N += 1
    return (suma/N)**0.5

def imprimir_resultados(resultados_rmsd):
    """
    Función para imprimir los resultados de RMSD y p-value en un archivo xlsx
    """
    import xlsxwriter
    # crear directorio resultados si no existe
    if not os.path.exists(os.path.abspath('.') + '/resultados'):
        os.makedirs(os.path.abspath('.') + '/resultados')
    xls_name = os.path.abspath('.') + '/resultados/resultados_rmsd.xlsx'
    workbook = xlsxwriter.Workbook(xls_name)
    worksheet = workbook.add_worksheet()
    row = 0
    worksheet.write(row, 0, "Metodo")
    worksheet.write(row, 1, "PDB")
    worksheet.write(row, 2, "RMSD")
    row = 1
    for metodo in resultados_rmsd:
        for pdb in resultados_rmsd[metodo]["pdb"]:
            worksheet.write(row, 0, metodo)
            worksheet.write(row, 1, pdb)
            worksheet.write(row, 2, resultados_rmsd[metodo]["pdb"][pdb]["rmsd"])
            row += 1
        worksheet.write(row, 1, "RMSD general")
        worksheet.write(row, 2, resultados_rmsd[metodo]["rmsd"])
        row += 2
    workbook.close()

def cargar_estructura(archivo_pdb):
    """
    Función para cargar estructuras desde archivos PDB
    """
    parser = PDBParser(QUIET=True)
    return parser.get_structure('estructura', archivo_pdb)

def get_atom_identifier(atom):
    """
    Función para obtener un identificador único para un átomo basado
    en el tipo de residuo, número de residuo y nombre del átomo.
    """
    residue_id = atom.get_parent().id[1]
    residue_name = atom.get_parent().resname
    atom_name = atom.get_id()
    return f"{residue_name}_{residue_id}_{atom_name}"

# Ruta de lista de archivos PDB a comparar
# deben tener la misma cantidad de pdbs y estar ordenados uno a uno
# se compara el pdb1 con el pdb1, pdb2 con pdb2, etc
# file_list_native_path = os.path.abspath('.')+"/pdb_list_origin_XRAY"
# file_list_processed_path = os.path.abspath('.')+"/pdb_list_processed_XRAY"
file_list_native_path = os.path.abspath('.')+"/pdb_list_origin_NMR"
file_list_processed_path = os.path.abspath('.')+"/pdb_list_processed_NMR_bp"


# inicializar variables
file_list_native = [] # lista de archivos que continen lista de pdbs
pdb_list_native = [] # lista de archivos pdb
file_list_processed = [] # lista de archivos que continen lista de pdbs
resultados_rmsd = {} # resultados de rmsd

# Leer lista de archivos PDB
with open(file_list_native_path) as f:
    file_list_native = (f.read().splitlines())
# Leer lista de archivos PDB del pdb_list
for path in file_list_native:
    with open(os.path.abspath('.')+"/"+path) as f:
        pdb_list_native = f.read().splitlines()

with open(file_list_processed_path) as f:
    file_list_processed = f.read().splitlines()

# recorremos la lista de resultados y leemos los archivos pdb uno a uno
for line in file_list_processed:
    metodo = line.split('/')[-1]
    print("Resultados RMSD para el archivo: ", metodo)

    # leemos los path pdb procesados
    pdb_list_processed = [] # lista de archivos pdb
    with open(os.path.abspath('.')+"/"+line) as f:
        pdb_list_processed = f.read().splitlines()
    if len(pdb_list_native) != len(pdb_list_processed):
        print("Error: las listas de pdbs no tienen la misma cantidad de pdbs")
        exit()

    atoms_structuras_nativas = []  # lista de atomos de todas las estructuras nativas
    atoms_structuras_procesadas = []  # lista de atomos de todas las estructuras procesadas ya sobrepuestas

    for i in range(len(pdb_list_native)):
        # ignoramos pdbs que no se pudieron procesar
        if pdb_list_processed[i] == '***':
            continue
        # Rutas de los archivos PDB de las dos estructuras a comparar
        archivo_pdb1 = os.path.abspath('.') +"/"+ pdb_list_native[i]
        archivo_pdb2 = os.path.abspath('.') +"/"+ pdb_list_processed[i]

        # Cargar las estructuras desde los archivos PDB
        estructura1 = cargar_estructura(archivo_pdb1)
        estructura2 = cargar_estructura(archivo_pdb2)

        # Obtener listas de átomos para superponer
        atoms_structure1 = estructura1.get_atoms()
        atoms_structure2 = estructura2.get_atoms()

        # Crear diccionarios para almacenar los átomos únicos de cada estructura
        unique_atoms_structure1 = {}
        unique_atoms_structure2 = {}

        # Obtener identificadores únicos de átomos para la estructura 1
        for atom in atoms_structure1:
            atom_id = get_atom_identifier(atom)
            unique_atoms_structure1[atom_id] = atom

        # Obtener identificadores únicos de átomos para la estructura 2
        for atom in atoms_structure2:
            atom_id = get_atom_identifier(atom)
            unique_atoms_structure2[atom_id] = atom

        # Crear listas para almacenar los átomos equivalentes entre ambas estructuras
        equivalent_atoms_structure1 = []
        equivalent_atoms_structure2 = []

        # Encontrar los átomos comunes entre ambas estructuras y almacenarlos en las listas respectivas
        for atom_id, atom in unique_atoms_structure1.items():
            if atom_id in unique_atoms_structure2:
                equivalent_atoms_structure1.append(atom)
                equivalent_atoms_structure2.append(unique_atoms_structure2[atom_id])

        # Inicializar el superimposer
        superimposer = Superimposer()

        # Establecer los átomos a superponer
        # aqui se sobre pone (rotar y transladar las estructuras)
        superimposer.set_atoms(equivalent_atoms_structure1, equivalent_atoms_structure2)

        # Realizar la superposición
        # aqui se aplica y se modifican las cordenadas
        superimposer.apply(equivalent_atoms_structure2)

        # cargar equivalent_atoms_structure en una sola estructura para 1 y 2
        atoms_structuras_nativas.append(equivalent_atoms_structure1)
        atoms_structuras_procesadas.append(equivalent_atoms_structure2)

        # Calcular el RMSD
        #rmsd = superimposer.rms
        rmsd = calcular_rmsd_ecuacion([equivalent_atoms_structure1], [equivalent_atoms_structure2])
        if estructura1.header['idcode']:
            pdb_name = estructura1.header['idcode']
        else:
            # scar de pdb_list_native[i]
            pdb_name = pdb_list_native[i].split('/')[-1].split('.')[0]
        print(f"Comparando pdb = {pdb_name}")
        print(f"El RMSD entre las dos estructuras es: {rmsd:.8f} Å")

        # Almacenar los resultados
        if metodo not in resultados_rmsd:
            # se propone guardar en una estrucutra definida por metodo de comparacion en el que se almacenara una lista
            # de pdbs con su rmsd y p_value ademas de tener un rmsd general de todos los pdbs
            resultados_rmsd[metodo] = {"pdb": {}, "rmsd": 0}
        if pdb_name not in resultados_rmsd[metodo]["pdb"]:
            resultados_rmsd[metodo]["pdb"][pdb_name] = {"rmsd": rmsd}
        else:
            print("Error: el pdb ya fue procesado. PDB duplicado")

    # Calcular el RMSD general entre todas las estructuras acumuladas
    rmsd_general = calcular_rmsd_ecuacion(atoms_structuras_nativas, atoms_structuras_procesadas)
    # almacenar el rmsd general
    resultados_rmsd[metodo]["rmsd"] = rmsd_general

# imprimir resultados en un archivo xlsx
imprimir_resultados(resultados_rmsd)



