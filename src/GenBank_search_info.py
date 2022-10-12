'''
NAME
    GenBank_search_info.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
    Programa que obtiene informacion de un archivo de GenBank y busca una lista 
    de genes dada por el usuario para regresar sus primeros 15 nucleotidos de 
    ADN, ARN y la secuencia proteica correspondiente

CATEGORY 
    GenBank
    
USAGE
    py .\src\GenBank_search_info.py -f path/to/file -o OUTPUT -g GENES 
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -f FILE, --File FILE
                        GenBank file
    -o OUTPUT, --output OUTPUT
                        Output file
    -g GENES, --genes GENES
                        List of genes
                       
INPUT
    Archivo GenBank
    
SEE ALSO
    GenBank_info.py
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
                    required=True)
                    
parser.add_argument("-g", "--genes",
                    help="Lista de genes a buscar",
                    type = str,
                    nargs = '+',
                    required=True)
                    
args = parser.parse_args()

# Declarar la funcion a usar
def GenBank_info(output, input, genes):
    
    # Abrir el archivo output
    with open(output, "w") as file:
        
        # Ciclo que parsea el archivo GenBank
        for gb_record in SeqIO.parse(input, "genbank"):
            
            # Escribir en el archivo output la informacion general del archivo
            file.write("DATABANK FILE INFO:\n\n")
            file.write(f"Organism:\n\t{gb_record.annotations['organism']}\n\n")
            file.write(f"Date:\n\t{gb_record.annotations['date']}\n\n")
            file.write(f"Country:\n\t{gb_record.features[0].qualifiers['country']}\n\n")
            file.write(f"Isolation Source:\n\t{gb_record.features[0].qualifiers['isolation_source']}\n\n")
            file.write("\nGENES INFO:\n\n")
            
            # Buscar cada gen de la lista en el archivo GenBank
            for gen in genes:
                for feature in gb_record.features:
                    if feature.type == "gene":
                        if feature.qualifiers['gene'][0] == gen:
                            
                            # Obtener los datos del gen requeridos
                            dna = gb_record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
                            dna = dna[:15]
                            rna = dna[:15].transcribe()
                            protein = dna[:15].translate()
                            
                            # Escribir los datos en el archivo output
                            file.write(f"Gen:\n\t{gen}\n\n")
                            file.write(f"\tDNA:\n\t\t{dna}\n\n")
                            file.write(f"\tRNA:\n\t\t{rna}\n\n")
                            file.write(f"\tProtein:\n\t\t{protein}\n\n")
    return(f"\nSe ha generado el archivo {output} con la informacion obtenida\n")
    
# Llamar a la funcion para obtener el archivo output
result_file = GenBank_info(args.output, args.file, args.genes)

# Indicar al usuario que se genero el archivo
print(result_file)

                        
                        