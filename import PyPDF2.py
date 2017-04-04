import PyPDF2

# file = open('experimento_1.pdf', 'rb')
#file = open('experimento_2.pdf', 'rb')
#file = open('sample.pdf', 'rb')
#file = open('Curso_Termodinamica_Python.pdf', 'rb')

file = open('Modeling.pdf', 'rb')
new_file = open('texto.txt', 'w')

file_pdf = PyPDF2.PdfFileReader(file)
number_pages = file_pdf.getNumPages()
pageObj = file_pdf.getPage(1)
information = pageObj.extractText()

#print(number_pages)
#print(pageObj)
#print(information[5302:])
print(information)

palabra = "Table1"

condicion_1 = palabra in information
print(condicion_1)

new_file.write(information)

file.close()
new_file.close()




