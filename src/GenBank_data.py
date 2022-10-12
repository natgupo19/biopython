'''
NAME
    GenBank_data.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
   Programa que obtiene informacion de los features y las anotaciones 
   de un archivo GenBank.

CATEGORY 
    GenBank
    
USAGE
    py .\src\GenBank_data.py -f path/to/file -o path/to/file 
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -f FILE, --file FILE
                        GenBank file
    -o OUTPUT, --output OUTPUT
                        Output file
                       
INPUT
    Archivo GenBank
    
SEE ALSO
    quality_records
'''

# Importar librerias
from Bio import SeqIO
import argparse

# Agregar el parser
parser = argparse.ArgumentParser(description = "Obtener datos de anotacion y features")

parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo GenBank",
                    required=True)
                                    
parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=False)
                    
args = parser.parse_args()

# Ciclo que parsea el archivo GenBank
for gb_data in SeqIO.parse(args.file, "genbank"):
    
    # Si el usuario desea obtener un archivo output...
    if args.output:
        
        # Abrir el archivo y escribir la informacion de las anotaciones y los features
        with open(args.output, "w") as file:
            file.write("Annotations:\n\n")
            file.write(f"Name:\n\t{gb_data.name}\n\n")
            file.write(f"Date:\n\t{gb_data.annotations['date']}\n\n")
            file.write(f"Organism:\n\t{gb_data.annotations['organism']}\n\n")
            file.write(f"\nFeatures:\n\n")
            file.write(f"Country:\n\t{gb_data.features[0].qualifiers['country']}\n\n")
            file.write(f"Molecular Type:\n\t{gb_data.features[0].qualifiers['mol_type']}\n\n")
            file.write(f"Isolation Source:\n\t{gb_data.features[0].qualifiers['isolation_source']}\n\n")
        
        # Indicar el path del archivo de salida
        print(f"\nLos datos se encuentran en el archivo: {args.output}\n")
    
    # Si no se pide archivo de salida...       
    else:
        
        # Imprimir la informacion de las anotaciones y los features
        print("\nAnnotations:\n")
        print(f"Name:\n\t{gb_data.name}\n")
        print(f"Date:\n\t{gb_data.annotations['date']}\n")
        print(f"Organism:\n\t{gb_data.annotations['organism']}\n")
        print(f"\nFeatures:\n")
        print(f"Country:\n\t{gb_data.features[0].qualifiers['country']}\n")
        print(f"Molecular Type:\n\t{gb_data.features[0].qualifiers['mol_type']}\n")
        print(f"Isolation Source:\n\t{gb_data.features[0].qualifiers['isolation_source']}\n")