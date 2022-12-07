'''
NAME
    Entrez_efetch_elink.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
    Programa que genera un termino para ser buscado por esearch y egquery, y así generar un archivo output con bases de datos asociados a sus ID's

CATEGORY 
    BioPython
    
USAGE
    py .\src\Entrez_efetch_elink.py -t TERMINO -o OUTPUT  
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -t TERM, --File File
                        Archivo que contiene los ID's a buscar
    -o OUTPUT, --output OUTPUT
                        Output file
                       
SEE ALSO
    Entrez_search_egquery.py
'''

from Bio import Entrez
import argparse
import re

# Agregar el parser
parser = argparse.ArgumentParser(description = "Obtener datos de anotacion y features")

parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo con ID's a buscar",
                    required=True)
                                    
parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=True)

args = parser.parse_args()

Entrez.email = "natgupo@lcg.unam.mx"

def paper_names(input, output):
    '''
    '''
    with open(input, 'r') as file:
        info = file.read()
    my_list = info.split("*")
    my_list = [element.split("\n") for element in my_list]
    my_list.pop(0)
    delete_space = [element.pop(-1) for element in my_list]
    
    with open(output, "w") as output:
        for data in my_list:
            output.write(f"Organismo: {data[0]}\n")
        
            Db_IDs = []
            for data_base in data[1:]:
                required_db = data_base.split(":")
                if required_db[0] == "pubmed" or required_db[0] == "pmc":
                    Db_IDs.append(required_db)
                    for id in Db_IDs[1:]:
                        handle = Entrez.efetch(db = required_db[0], id = id, rettype = "medline", retmode = "text")
                        title = handle.read()
                        citas = Entrez.read(Entrez.elink(dbfrom = required_db[0], db = required_db[0], LinkName = "pubmed_pmc_refs", from_uid = id))
                        handle.close()
                        lines = str(title).split("\n")
                        for line in lines:
                            i = re.search("TI", line)
                            if i:
                                output.write("Titulo: " + str(i) + "\n\n")
                                
                  
                        
                
    
    return()
    
miau = paper_names(args.file, args.output)
print(miau)
    
        
    

