#!/usr/bin/env python

# Copyright 2023, Hector Diaz
# This program is distributed under General Public License v. 3.
# COPYING for a copy of the license.

__description__ = \
    """
    Determines the absolute precision between residues.
    """

__author__ = "Hector Diaz"
__date__ = "08092021"

import os, configparser
from lib.funtions import print_to_excel, search_line
from lib.pdb_torsionCHI import torsionCHI


def main():
    """
    Call if this is called from the command line.
    Calcula el porcentaje de residuos correctos segun la medida de precision absoluta basada en un angulo de aceptacion
    ejecucion: python3 main.py
    """
    # leer archivos con lista de ubicaciones de pdbs el primero debe contener a los originales
    # se espera que las listas esten ordenadas y tengan la misma cantidad de pdbs, en el caso que un pdb sea procesable
    # se debe colocar *** en la lista en el caso en que un pdb sea improcesable
    # todos los pdbs deben estar previamente estructurados de manera a tener la misma secuencia de residuos
    print('#######################################################')
    print("///// Resultados de medida por Exactitud Absoluta /////")
    print('#######################################################')
    res_metodos = {}  # diccionario para almacenar los resultados de cada metodo
    config = configparser.ConfigParser()
    config.read(os.path.abspath('.')+"/config.ini")
    correct_angle = config.getint('configuracion', 'correct_angle')  # angulo tomado como correcto
    print_pdb_console = config.getboolean('configuracion', 'print_pdb_console')  # imprimir pdbs en consola
    print_pdb_excel = config.getboolean('configuracion', 'print_pdb_excel')  # imprimir pdbs en excel
    pdb_origin_file = os.path.abspath('.')+"/"+config.get('configuracion','pdb_origin_file') # archivo con lista de pdbs originales
    pdb_processed_files = os.path.abspath('.')+"/"+config.get('configuracion','pdb_processed_files') # archivo con lista de pdbs procesados

    file_list = []
    with open(pdb_origin_file, 'r') as list_origin:
        file_list.append(list_origin.read().splitlines())
    with open(pdb_processed_files, 'r') as list_processed:
        file_list.append(list_processed.read().splitlines())
    pdbs_file_list = {}  # lista de listas a ser procesados
    pdbs_file_key = []  # nombre de cada lista
    for l in range(len(file_list)):
        for pdb_list_file in file_list[l]:
            list_name = pdb_list_file.split('/')[-1]
            pdbs_file_key.append(list_name)
            f = open(os.path.abspath('.')+"/"+pdb_list_file, 'r')
            pdbs_file_list[list_name] = f.read().splitlines()
            f.close()

    # recorrer por listas para comparar la primera con cada una de las restantes
    # k es el identificador del archivo lista a comparar parte desde 1 por que el 0 esta reservado para la lista origen
    for k in range(1, len(pdbs_file_key)):
        # recorrer dentro de las listas para obtener los pdbs a comparar
        # h es el identificador del pdb dentro del archivo lista a comparar
        for h in range(len(pdbs_file_list[pdbs_file_key[k]])):
            # por cada lista se tiene un pdb
            pdb_file = pdbs_file_list[pdbs_file_key[0]][h]
            if pdb_file[:3] == "***":
                continue
            # almacenamos los rotameros del pdb original
            rotameros_nativo = torsionCHI(pdb_file)
            pdb_file = pdbs_file_list[pdbs_file_key[k]][h]
            # se agrega un filtro para saltar la comparacion si se encuentra una linea de *** en la lista de pdb
            if pdb_file[:3] == "***":
                continue
            # almacenamos los rotameros del pdb a comparar
            rotameros_metodo = torsionCHI(pdb_file)

            correctox1 = 0  # cantidad de correctos X1
            evaluadox1 = 0  # cantidad de residuos validos para evaluar x1
            correctox1x2 = 0  # cantidad de correctos X1X2
            evaluadox1x2 = 0  # cantidad de residuos validos para evaluar x1x2

            # agregar columnas
            rotameros_metodo[0] = rotameros_metodo[0][:-1] + "%10s%10s%15s%15s%10s%10s" % (
                "Ch1orig", "Ch2orig", "EvaluableX1", "EvaluableX1+2", "X1", "X1+2\n")

            # out.append("%30s%4s \"%s\"%10.2F%10.2F%10.2F%10.2F%10.2F%10.2F%10.2F\n" % (short_pdb, labels[i][:3], labels[i][4:], phi_psi[i][0], phi_psi[i][1],dihedrals[i][0], dihedrals[i][1], dihedrals[i][2], dihedrals[i][3],dihedrals[i][4]))

            for j in range(1, len(rotameros_metodo)):
                # almacenamos los chi1 y chi2 originales y predichos
                # siendo el pdb native (rotameros_nativo) el pdb de origen
                # buscamos la linea que pertenece al mismo residuo en el pdb native para obtener los valores de chi1 y chi2
                linea_nativo = search_line(rotameros_nativo, rotameros_metodo[j][47:51], rotameros_metodo[j][46:47])
                if not linea_nativo:
                    print
                    "Error: no se encontro el rotamero en el pdb nativo" + rotameros_metodo[j][10:51]
                    continue
                ori_amino = linea_nativo[41:44]
                ori_chi1 = float(linea_nativo[72:82])
                ori_chi2 = float(linea_nativo[82:92])
                pred_amino = rotameros_metodo[j][41:44]
                pred_chi1 = float(rotameros_metodo[j][72:82])
                pred_chi2 = float(rotameros_metodo[j][82:92])
                rotameros_metodo[j] = rotameros_metodo[j][:-1] + "%10.2F%10.2F" % (ori_chi1, ori_chi2)

                # Calculo de diferencia entre angulos teniendo en cuenta el cuadrante en el que se encuentra.
                Sx1 = 360 # valor de diferencia de angulo para x1, se inicializa en 360 para que sea mayor a cualquier valor
                Sx1x2 = 360 # valor de diferencia de angulo para x1+2, se inicializa en 360 para que sea mayor a cualquier valor
                # consideraciones para evaluar x1 y x1x2
                if ori_chi1 != 0 and pred_chi1 != 0:
                    # Se resta el valor menor al mayor para obtener el valor absoluto de la diferencia
                    Ex1 = max(ori_chi1, pred_chi1) - min(ori_chi1, pred_chi1)
                    # Se evalua si el valor absoluto de la diferencia es mayor a 180
                    Sx1 = min(Ex1, 360 - Ex1)
                    evaluadox1 = evaluadox1 + 1  # aumentar la cantidad de residuos evaluados
                    rotameros_metodo[j] = rotameros_metodo[j] + "%15s" % ("1")  # agregar el valor para la tabla
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
                            Sx1x2 = min(dif1, dif2)
                        else:
                            Ex1x2 = max(ori_chi2, pred_chi2) - min(ori_chi2, pred_chi2)
                            Sx1x2 = min(Ex1x2, 360 - Ex1x2)
                        evaluadox1x2 = evaluadox1x2 + 1  # aumentar la cantidad de residuos evaluados
                        rotameros_metodo[j] = rotameros_metodo[j] + "%15s" % ("1")  # agregar el valor para la tabla
                    else:
                        rotameros_metodo[j] = rotameros_metodo[j] + "%15s" % ("0")  # agregar el valor para la tabla
                else:
                    rotameros_metodo[j] = rotameros_metodo[j] + "%15s%15s" % (
                        "0", "0")  # agregar el valor para la tabla
                if Sx1 <= correct_angle:
                    rotameros_metodo[j] = rotameros_metodo[j] + "%10s" % ("1")
                    correctox1 = correctox1 + 1  # aumentar la cantidad de correctos para el calculo de porcentaje
                    if Sx1x2 <= correct_angle:
                        rotameros_metodo[j] = rotameros_metodo[j] + "%10s" % ("1\n")
                        correctox1x2 = correctox1x2 + 1
                    else:
                        rotameros_metodo[j] = rotameros_metodo[j] + "%10s" % ("0\n")
                else:
                    rotameros_metodo[j] = rotameros_metodo[j] + "%10s%10s" % ("0", "0\n")

            prob_x1 = correctox1 * 100 / evaluadox1 # probabilidad de correctos para x1
            prob_x1x2 = correctox1x2 * 100 / evaluadox1x2 # probabilidad de correctos para x1+2
            # impresion de resultados
            estado_arte = pdbs_file_key[k]
            print("Estado del arte: %s" % (estado_arte))
            if pdbs_file_list[pdbs_file_key[0]][h].split('/'):
                if pdbs_file_list[pdbs_file_key[0]][h].split('/')[-1].split(".")[0]:
                    pdb_name = pdbs_file_list[pdbs_file_key[0]][h].split('/')[-1].split(".")[0]
                else:
                    pdb_name = pdbs_file_list[pdbs_file_key[0]][h].split('/')[-1]
            else:
                pdb_name = pdbs_file_list[pdbs_file_key[0]][h]
            print("PDB: %s" % (pdb_name))
            print("%20s%20s" % ("X1(%)", "X1+2(%)"))
            print("%20.2F%20.2F" % (prob_x1, prob_x1x2))
            print("%10s%10.0F%20.0F" % ("Evaluado  ", evaluadox1, evaluadox1x2))
            print("%10s%10.0F%20.0F" % ("Correcto  ", correctox1, correctox1x2))

            # Sumatorias por estado del arte
            if not estado_arte in res_metodos:
                res_metodos[estado_arte] = {}
                res_metodos[estado_arte] = {
                    'metodo': estado_arte,
                    'residuos': len(rotameros_metodo) - 1,
                    'evaluadox1': evaluadox1,
                    'evaluadox1x2': evaluadox1x2,
                    'correctox1': correctox1,
                    'correctox1x2': correctox1x2,
                }
            else:
                res_metodos[estado_arte] = {
                    'metodo': estado_arte,
                    'residuos': res_metodos[estado_arte]['residuos'] + len(rotameros_metodo) - 1,
                    'evaluadox1': res_metodos[estado_arte]['evaluadox1'] + evaluadox1,
                    'evaluadox1x2': res_metodos[estado_arte]['evaluadox1x2'] + evaluadox1x2,
                    'correctox1': res_metodos[estado_arte]['correctox1'] + correctox1,
                    'correctox1x2': res_metodos[estado_arte]['correctox1x2'] + correctox1x2,
                }
            # imprimir en consola
            if print_pdb_console:
                print()
                print("".join(rotameros_metodo))
            # imprimir en excel los resultados
            if print_pdb_excel:
                print_to_excel(rotameros_metodo, pdb_name, estado_arte)

    # imprimir res_metodos
    for est_arte in res_metodos:
        print('Metodo: %s' % (res_metodos[est_arte]['metodo']))
        print('%15s%15s%15s%15s%25s' % ('Residuos', 'Evaluado X1', 'Evaluado X1+2', 'Correcto X1', 'Correcto X1+2'))
        print('%15.0F%15.0F%15.0F%15.0F%15.0F' % (
            res_metodos[est_arte]['residuos'], res_metodos[est_arte]['evaluadox1'],
            res_metodos[est_arte]['evaluadox1x2'],
            res_metodos[est_arte]['correctox1'], res_metodos[est_arte]['correctox1x2']))
        p1 = float(res_metodos[est_arte]['correctox1']) * 100 / float(res_metodos[est_arte]['evaluadox1'])
        p2 = float(res_metodos[est_arte]['correctox1x2']) * 100 / float(res_metodos[est_arte]['evaluadox1x2'])

        print('%15s%15.2F%17s%15.2F' % ('Probabilidad X1', p1, 'Probabilidad X1+2', p2))
        print('#############################################')
        print()


if __name__ == "__main__":
    main()