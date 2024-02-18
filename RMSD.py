import os
from Bio.PDB import PDBParser, Superimposer
from scipy.stats import ttest_rel


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
    worksheet.write(row, 3, "P-value")
    row = 1
    for metodo in resultados_rmsd:
        for pdb in resultados_rmsd[metodo]["pdb"]:
            worksheet.write(row, 0, metodo)
            worksheet.write(row, 1, pdb)
            worksheet.write(row, 2, resultados_rmsd[metodo]["pdb"][pdb]["rmsd"])
            # worksheet.write(row, 3, ', '.join([f'{p:.8f}' for p in resultados_rmsd[metodo]["pdb"][pdb]["p_value"]]))
            for j in range(len(p_value)):
                worksheet.write(row, 3+j, resultados_rmsd[metodo]["pdb"][pdb]["p_value"][j])
            row += 1
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
file_list_native_path = os.path.abspath('.')+"/pdb_list_origin_XRAY"
file_list_processed_path = os.path.abspath('.')+"/pdb_list_processed_XRAY"

# inicializar listas
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
        superimposer.set_atoms(equivalent_atoms_structure1, equivalent_atoms_structure2)

        # Realizar la superposición
        superimposer.apply(estructura2.get_atoms())

        # Calcular el RMSD
        rmsd = superimposer.rms
        print(f"Comparando pdb = {estructura1.header['idcode']}")
        print(f"El RMSD entre las dos estructuras es: {rmsd:.8f} Å")

        # Obtener las coordenadas de los átomos después de la superposición
        coords1 = [atom.get_coord() for atom in equivalent_atoms_structure1]
        coords2 = [atom.get_coord() for atom in equivalent_atoms_structure2]

        # Realizar la prueba t pareada para comparación de residuos
        t_statistic, p_value = ttest_rel(coords1, coords2)

        # print(f"Valor de significancia (P) para la comparación de residuos: {p_value:.4f}" .join)
        print(f"Valor de significancia (P) para la comparación de residuos: {', '.join([f'{p:.8f}' for p in p_value])}")

        # Almacenar los resultados
        if metodo not in resultados_rmsd:
            # se propone guardar en una estrucutra definida por metodo de comparacion en el que se almacenara una lista
            # de pdbs con su rmsd y p_value ademas de tener un rmsd y p_value general de todos los pdbs
            resultados_rmsd[metodo] = {"pdb": {}, "rmsd": 0, "p_value": []}
        if estructura1.header['idcode'] not in resultados_rmsd[metodo]["pdb"]:
            resultados_rmsd[metodo]["pdb"][estructura1.header['idcode']] = {"rmsd": rmsd, "p_value": p_value}
        else:
            print("Error: el pdb ya fue procesado. PDB duplicado")

imprimir_resultados(resultados_rmsd)


