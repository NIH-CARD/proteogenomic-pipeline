from PyPDF2 import PdfFileMerger
import os
import sys

#Create and instance of PdfFileMerger() class
merger = PdfFileMerger()

#Create a list with PDF file names
path_to_files = sys.argv[1] #r'pdf_files/'

#Get the file names in the directory
for root, dirs, file_names in os.walk(path_to_files):
    #Iterate over the list of file names
    for file_name in file_names:
        #Append PDF files
        merger.append(path_to_files + file_name)

#Write out the merged PDF
merger.write(sys.argv[1]+"all_sashimis.pdf")
merger.close()
