'''
NAME
    quality_records.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
   Programa que selecciona las secuencias en las que todos sus nucleotidos 
   superen el umbral de calidad para el valor de Qscore.

CATEGORY 
    DNA sequence
    
USAGE
    py .\src\quality_records.py -f path/to/file -o path/to/file -u UMBRAL
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -f FASTQ, --Fastq FASTQ
                        Fastq file with DNA sequence
    -o OUTPUT, --output OUTPUT
                        Output file
    -u UMBRAL, --umbral UMBRAL
                        Umbral value
                       
INPUT
    Archivo con secuencias de DNA a analizar
    
SEE ALSO
    max_orf_to_protein
    get_codons
'''

# Importar librerias
from Bio import SeqIO
import argparse

# Agregar el parser
parser = argparse.ArgumentParser(description = "Obtener secuencias que superen el umbral de calidad")

parser.add_argument("-f", "--Fastq",
                    metavar="path/to/file",
                    help="Archivo con secuencia de ADN",
                    required=True)
                                    
parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=True)
                    
parser.add_argument("-u", "--umbral",
                    type= int,
                    help="Valor de umbral minimo esperado",
                    required=True)
                    
args = parser.parse_args()

# Declarar la funcion que obtendra los datos
def umbral_quality(archivo, umbral):
    '''
    Seleccionar aquellas secuencias en las que todos sus nucleotidos 
    superen el umbral de calidad.
        Parameters:
            archivo (str): Path del archivo fastq
            umbral (int): Numero que determina el umbral de calidad 
        Returns:
            lectures (int): Numero de secuecias que pasan el umbral 
    '''
    # Abrir el archivo de salida
    with open(args.output, "w") as file:
        
        # Inicializar el contador
        lectures = 0
        
        # Ciclo para parsear la secuencia
        for record in SeqIO.parse(archivo, "fastq"):
            
            # Solo tomar en cuenta las secuencias que su qscore sea mayor o igual al umbral
            scores = record.letter_annotations["phred_quality"]
            if min(scores) >= umbral:
                
                # Escribir en el archivo el id y la secuencia
                file.write(f"{record.id}\n{record.seq}\n\n")
                
                # Contar las secuencias que pasen el umbral
                lectures += 1
    return(f"Numero de lecturas que pasan el umbral: {lectures}\n")

# Llamar a la funcion 
total_records = umbral_quality(args.Fastq, args.umbral)

# Imprimir el numero de lecturas encontradas y el path del archivo output
print(f"\n{total_records}")
print(f"Los datos se encuentran en el archivo: {args.output}\n")