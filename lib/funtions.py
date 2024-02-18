
import xlsxwriter
import os

# imprimir excel desde una lista de str formateando por celdas
def print_to_excel(rotameros, pdb_name, estado_arte, correct_angel):

    # crear directorio estado_arte si no existe
    if not os.path.exists(os.path.abspath('.')+'/resultados/'+estado_arte+'_'+str(correct_angel)):
        os.makedirs(os.path.abspath('.')+'/resultados/'+estado_arte+'_'+str(correct_angel))
    xls_name= os.path.abspath('.')+'/resultados/'+estado_arte+'_'+str(correct_angel)+'/pa_'+pdb_name+"_"+estado_arte+'.xlsx'
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

def init_excel_results():
    if not os.path.exists(os.path.abspath('.') + '/resultados'):
        os.makedirs(os.path.abspath('.') + '/resultados')
    xls_name = os.path.abspath('.') + '/resultados/resultados_precision_absoluta.xlsx'
    workbook = xlsxwriter.Workbook(xls_name)
    worksheet = workbook.add_worksheet()
    row = 0
    worksheet.write(row, 0, "Metodo")
    worksheet.write(row, 1, "PDB")
    worksheet.write(row, 2, "Residuos")
    worksheet.write(row, 3, "Evaluado X1")
    worksheet.write(row, 4, "Evaluado X1+2")
    worksheet.write(row, 5, "Correcto X1")
    worksheet.write(row, 6, "Correcto X1+2")
    worksheet.write(row, 7, "Probabilidad X1")
    worksheet.write(row, 8, "Probabilidad X1+2")
    return workbook, worksheet

def xls_insert_resul(worksheet, row, dict):
    worksheet.write(row, 0, dict['metodo'])
    worksheet.write(row, 1, dict['pdb'])
    worksheet.write(row, 2, dict['residuos'])
    worksheet.write(row, 3, dict['evaluadox1'])
    worksheet.write(row, 4, dict['evaluadox1x2'])
    worksheet.write(row, 5, dict['correctox1'])
    worksheet.write(row, 6, dict['correctox1x2'])
    p1 = float(dict['correctox1']) * 100 / float(dict['evaluadox1'])
    p2 = float(dict['correctox1x2']) * 100 / float(dict['evaluadox1x2'])
    worksheet.write(row, 7, p1)
    worksheet.write(row, 8, p2)


    return worksheet, row

def print_to_excel_results(res_metodos):
    # crear directorio resultados si no existe
    if not os.path.exists(os.path.abspath('.') + '/resultados'):
        os.makedirs(os.path.abspath('.') + '/resultados')
    xls_name = os.path.abspath('.') + '/resultados/resultados_precision_absoluta.xlsx'
    workbook = xlsxwriter.Workbook(xls_name)
    worksheet = workbook.add_worksheet()
    row = 0
    worksheet.write(row, 0, "Metodo")
    worksheet.write(row, 1, "Residuos")
    worksheet.write(row, 2, "Evaluado X1")
    worksheet.write(row, 3, "Evaluado X1+2")
    worksheet.write(row, 4, "Correcto X1")
    worksheet.write(row, 5, "Correcto X1+2")
    worksheet.write(row, 6, "Probabilidad X1")
    worksheet.write(row, 7, "Probabilidad X1+2")
    row = 1
    for est_arte in res_metodos:
        worksheet.write(row, 0, res_metodos[est_arte]['metodo'])
        worksheet.write(row, 1, res_metodos[est_arte]['residuos'])
        worksheet.write(row, 2, res_metodos[est_arte]['evaluadox1'])
        worksheet.write(row, 3, res_metodos[est_arte]['evaluadox1x2'])
        worksheet.write(row, 4, res_metodos[est_arte]['correctox1'])
        worksheet.write(row, 5, res_metodos[est_arte]['correctox1x2'])
        p1 = float(res_metodos[est_arte]['correctox1']) * 100 / float(res_metodos[est_arte]['evaluadox1'])
        p2 = float(res_metodos[est_arte]['correctox1x2']) * 100 / float(res_metodos[est_arte]['evaluadox1x2'])
        worksheet.write(row, 6, p1)
        worksheet.write(row, 7, p2)
        row += 1
    workbook.close()


def search_line(rotamero_nativo, residuo, cadena):
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
