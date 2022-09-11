'''
NAME
    ARN_to_protein
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
    Programa que transforma una secuencia de RNA a 
    una cadena lipidica

CATEGORY 
    RNA sequence and protein 
    
USAGE
    py .\src\ARN_to_protein.py [-f path/to/file] [-s] [-o]
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -f FILE, --file FILE
                        File with DNA sequence
    -s SEQUENCE, --sequence SEQUENCE
                        Minimum size to search
    -o OUTPUT, --output OUTPUT
                        Output file
                        
INPUT
    Secuencia de RNA a analizar
    
SEE ALSO
    None
'''

import re
import argparse

# Agregar el parser
parser = argparse.ArgumentParser(description = "")

parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo con secuencia de ARN",
                    required=False)
                                    
parser.add_argument("-s", "--sequence",
                    help="Secuencia de ARN",
                    type=str,
                    required=False)
                    
parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=False)
                    
args = parser.parse_args()

# Si el usuario ingresa un archivo, leemos la secuencia
if args.file:
    with open(args.file, 'r') as sequence:
        rna = sequence.read().rstrip('\n').upper()

# Si el usuario inserta una secuencia de RNA
if args.sequence:
    rna = args.sequence.upper()
    
# Definimos la funcion para validar la secuencia de RNA                        
def validate(rna):
    '''
    Evalua si el archivo contiene algun caracter que no sea un nucleotido
        Parameters:
            rna (str): secuencia de RNA a analizar
        Returns:
            dna (Any): secuencia de DNA a transformar
            1 (int): si la secuencia es correcta
    '''
    # Buscar los caracteres invalidos y contarlos
    invalid = re.finditer("[^AUCG]", rna)
    match = len([*invalid]) 
    if not match:
        dna = rna.replace("U","T")
        return(dna)
    # Notificar el error si es que se encuentran caracteres invalidos
    if match:
        print (f"\nNo se introdujo una secuencia de RNA\n")
        return(0)

# Definimos la funcion para convertir una secuencia de DNA a proteina      
def dna_to_protein(dna):
    '''
    Convierte una cadena de DNA a una cadena peptidica
        Parameters:
            dna (Any): secuencia de RNA a analizar
        Returns:
            aminoacids (str): secuencia de DNA a transformar
    '''
    # Definimos el diccionario de codones
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T',
    'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N',
    'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R',
    'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H',
    'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',
    'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V',
    'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G',
    'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S',
    'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L',
    'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    # Declaramos la cadena de aminoacidos
    aminoacids = "" 
    # Recorremos la cadena de DNA e identificamos los codones
    if dna:
        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]
            aminoacids += gencode.get(codon, "*(codon incompleto)")
    return(aminoacids)

# Obtenemos la cadena de DNA y la cadena proteica   
DNA = validate(rna)
protein = dna_to_protein(DNA)

# Si el usuario lo indica, obtener un archivo con el resultado
if args.output:
    with open(args.output, 'w') as output_file:
        print(f"\nSecuencia proteica encontrada es: {protein}\n",
                file=output_file)
    print(f"\nSe ha generado el archivo {args.output} con la secuencia proteica obtenida\n")

# Imprimir los resultados   
else:
    print(f"\nLa secuencia proteica obtenida es: {protein}\n")
