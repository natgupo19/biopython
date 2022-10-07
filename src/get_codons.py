'''
NAME
    get_codons.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
    Imprime los codones, separados por un espacio, para cada uno de los 
    6 marcos de lecturas de cada secuencia introducida en formato FASTA

CATEGORY 
    DNA sequence and codons 
    
USAGE
    py .\src\get_codons.py -f path/to/file -o output
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -f FILE, --file FILE
                        File with DNA sequence
    -o OUTPUT, --output OUTPUT
                        Output file
                       
INPUT
    Archivo con secuencias de DNA a analizar
    
SEE ALSO
    max_orf_to_protein
'''

# Importar librerias
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
import argparse
import re

# Agregar el parser
parser = argparse.ArgumentParser(description = "Obtener secuencia proteica con mas aminoacidos")

parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo con secuencia de ADN",
                    required=True)
                                    
parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=True)
                    
args = parser.parse_args()

# Parsear el archivo fasta en un diccionario
dictionary = SeqIO.to_dict(SeqIO.parse(args.file, "fasta"))

# Abrir el archivo fasta
with open(args.output, "w") as file:
    
    # Buscar el primer marco de lectura
    for id in dictionary:
        
        # Buscar codon de inicio y checar que haya marco de lectura
        start_codon = SeqUtils.nt_search(str(dictionary[id].seq), "ATG")
        
        # Indicar al usuario si no se encuentra el marco de lectura
        if len(start_codon) == 1:
            print("No se encontraron orfs.")
        
        # Si no se encuentra el marco, seguir con el siguiente record
        else:
            
            # Declarar los valores y las condiciones necesarias para ir cambiando el marco de lectura
            i = 0
            while i <= 1:
                j = 1
                while j <= 3:
                    
                    # Obtener los codones de los primeros tres marcos de lectura
                    # Escribir los codones en el archivo
                    if(i == 0):
                        file.write(f">{id} en marco de lectura #{j}\n")
                        codons = []
                        for codon in re.findall(r"(.{3})", str(dictionary[id].seq[j-1::])):
                            codons.append(codon)
                        file.write(" ".join(codons))
                        file.write("\n\n")
                        
                    # Obtener los codones de los tres marcos de lectura de la secuencia complementaria
                    # Escribir los codones en el archivo
                    else:
                        file.write(f">{id} en marco de lectura # -{j}\n")
                        complements = []
                        for complement in re.findall(r"(.{3})", str(((dictionary[id].seq).complement())[j-1::])):
                            complements.append(complement)
                        file.write(" ".join(complements))
                        file.write("\n\n")
                    
                    # Aumentamos el valor de j e i para que cambie de marco de lectura
                    j += 1
                i += 1