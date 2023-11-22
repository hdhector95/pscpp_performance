
import xlsxwriter
import os

# imprimir excel desde una lista de str formateando por celdas
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