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
    Programa que busca orfs en una secuencia, obtiene sus secuencias proteicas 
    y selecciona la que tiene una mayor cantidad de aminoacidos

CATEGORY 
    DNA sequence and codons 
    
USAGE
    py .\src\get_codons.py [-f] [-o]
    
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

