'''
NAME
    Entrez_search_egquery.py
    
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
    py .\src\Entrez_search_egquery.py -t TERMINO -o OUTPUT  
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -t TERM, --term TERM
                        Termino a introducir
    -o OUTPUT, --output OUTPUT
                        Output file
                       
SEE ALSO
    Entrez_efetch_elink.py
'''

from Bio import Entrez
import argparse

# Agregamos el parser
parser = argparse.ArgumentParser(description="Obtener IDs de las bases de datos donde se encuentre la informacion dada")

parser.add_argument("-t", "--term",
                    help="Introduzca el organismo y los genes en el siguiente formato: 'Organismo1: Gene1, Gene2; Organismo2: Gene3, Gene4...'",
                    type=str,
                    required=True)

parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=True)

args = parser.parse_args()

Entrez.email = "natgupo@lcg.unam.mx"

# Definir la funcion para obtener el termino a buscar
def create_term(info):
    '''
    Funcion que toma un string con formato: "Organismo1: Gene1, Gene2; Organismo2: Gene3, Gene4..." para generar terminos que puedan utilizarse en egquery para buscar genes asociados a organismos.
        Parameters:
            info (str): La informacion del termino que se busca generar
        Returns:
            terms (list): Lista de los terminos generados
    '''
    # Inicializar la lista en donde se almacenaran los terminos
    terms = []
    
    # Separar la informacion mediante splits para generar una lista
    # con los datos deseados, en el formato requerido
    # Cada elemento de la lista correspondera a un organismo y sus respectivos genes
    # Inicializar un ciclo for para separar a cada organismo con sus genes
    for inf in info.split('; '):
        inf = inf.split(':')
        
        # Introducir el organismo en la lista
        term = inf[0] + "[Orgn] AND ("
        i = True
        
        # Separar los genes e introducirlos en el formato requerido
        for gen in inf[1].replace(' ', '').split(','):
            if i:
                term += gen + "[Gene]"
            else:
                term += " OR " + gen + "[Gene]"
            i = False
        term += ")"
        
        # Aniadir cada termino a la lista
        terms += [term] 
        
        
    return(terms)

# Definir la funcion para obtener los ID's
def term_ids(terms, path):
    '''
    '''
    # Abrir el archivo que contendra la informacion
    with open(path, "w") as file:
        for termino in terms:
            file.write("\n")
            organism = termino.split("[Orgn]")
            organism = organism[0]
            
            # Utilizar el termino en egquery
            handle = Entrez.egquery(term = termino)
            record = Entrez.read(handle)
            
            # Escribir el nombre del organismo en el archivo
            file.write(f"{organism}\n")
            
            # Inicializar el vector donde se guardaran las bases de datos
            dbases = []
            
            # Buscamos en los row las bases de datos que sí se encuentra
            for row in record["eGQueryResult"]:
                if row["Count"] != "0" and row["Count"] != "Error":
                    
                    # Aniadimos las bases de datos a una lista
                    dbases += [row["DbName"]]
            
            for dba in dbases:
                
                # Inicializar el vector donde se guardaran los ID's
                ids = []
                
                # Utilizar esearch para buscar los terminos en cada base de datos encontrada
                db_record = Entrez.read(Entrez.esearch(db = dba, term = termino))
                
                # Dar formato a la informacion 
                ids = db_record["IdList"]
                dba = dba.split(",")
                db_id = dba[0] + ": "
                db_ids = db_id + str(ids).replace('[', '').replace(']', '').replace("'", "")
                
                # Escribir la informacion en el archivo
                file.write(f"\n{db_ids}\n")
    return()

# Llamar a las funciones y obtener el archivo output
termino = create_term(args.term)
data_bases = term_ids(termino, args.output)

