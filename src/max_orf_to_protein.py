'''
NAME
    max_orf_to_protein.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
    Programa que busca orfs en una secuencia, obtiene sus secuencias proteicas 
    y selecciona la que tiene una mayor cantidad de aminoacidos

CATEGORY 
    DNA sequence and proteins 
    
USAGE
    py .\src\max_orf_to_protein.py [-s] [-f] [-o]
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -s SEQUENCE, --sequence SEQUENCE
                        DNA sequence by input
    -o OUTPUT, --output OUTPUT
                        Output file
    -f FILE, --file FILE
                        File with DNA sequence
                        
INPUT
    Secuencia de DNA a analizar
    
SEE ALSO
    ARN_to_protein
'''

# Importando librerias
import Bio
from Bio.Seq import Seq, MutableSeq
from Bio.SeqUtils import nt_search
import argparse
import warnings

# Agregar el parser
parser = argparse.ArgumentParser(description = "Obtener secuencia proteíca con mas aminoacidos")

parser.add_argument("-s", "--sequence",
                    help="Secuencia de ADN",
                    type=Seq,
                    required=False)

parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo con secuencia de ADN",
                    required=False)
                                    
parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=False)
                    
args = parser.parse_args()

# Ignorando warnings de Biopython
warnings.filterwarnings('ignore')
warnings.warn('BiopythonWarning') 

# Crear la funcion que obtendra la secuencia proteica mas larga
def largest_protein(dna): 
    '''
    Toma una secuencia de DNA y la traduce a una secuencia proteica a partir de los orfs 
    encontrados, regresa el peptido con mayor cantidad de aminoacidos. 
        Parameters:
            dna (Seq): Secuencia de nucleotidos
        Returns:
            protein (Seq): Secuencia proteica con mayor cantidad de aminoacidos
    '''   
    # Obtener las posiciones donde empiezan los orfs
    posiciones = nt_search(str(dna), "ATG")

    # Si no se encuentran orfs, indicar al usuario
    if len(posiciones) == 1:
        return("\nNo se encontraron orfs en la secuencia.\n")
    
    # Si se encuentran orfs...
    else:
        
        # Obtenemos todas las cadenas peptidicas por medio de los orfs
        peptides = []
        for i in range(1, len(posiciones)):
            start_codon = posiciones[i]
            fragments = dna[start_codon:]
            peptide = fragments.translate(to_stop = True)
            peptides.append(str(peptide))
        
        # Obtenemos la secuencia proteica mas larga
        lenghts = []    
        for size in peptides:
            lenghts.append(len(size))
            protein = peptides[lenghts.index(max(lenghts))]
        return(protein)

# Si la secuencia se introduce por input
if args.sequence:
    protein = largest_protein(args.sequence)
    print(f"\nLa secuencia proteica con mayor contenido de aminoacidos es: {protein}\n")

# Si la secuencia se encuentra en un archivo 
if args.file:
    with open(args.file, 'r') as archivo:
        ADN = archivo.read().replace('\n','')
    protein = largest_protein(Seq(ADN))
    print(f"\nLa secuencia proteica con mayor contenido de aminoacidos es: {protein}\n")
    
# Si se quiere pasar el resultado a un archivo
if args.output:
    with open(args.output, 'w') as archivo:
        archivo.write(f"La secuencia proteica con mayor contenido de aminoacidos es:\n{protein}")
    print(f"\nSe ha generado el archivo {args.output} con la secuencia proteica obtenida\n")

    
    




