'''
NAME
    AT_rich
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
     

CATEGORY
    DNA sequence
    
USAGE
    py .\src\AT_rich.py -f path/to/file [-s]
    
ARGUMENTS
    -h, --help          Show this help message and exit
    -f FILE, --file FILE
                        File with DNA sequence
    -s SIZE, --size SIZE
                        Minimum size to search
    -o path/to/output/file, --output path/to/output/file
                        Path for the output file
    -r ROUND, --round ROUND
                        Number of digits to round
                        
INPUT
    Secuencia de DNA a analizar
    
SEE ALSO
    None
'''

import argparse
import re 

parser = argparse.ArgumentParser(description = "Buscar regiones ricas en AT")

parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo con secuencia de ADN",
                    required=True)
         
parser.add_argument("-s", "--size",
                    help="cantidad minima de AT a buscar",
                    type=int,
                    required=False)

parser.add_argument("-o", "--output",
                    metavar="path/to/output/file",
                    help="Path for the output file",
                    required=False)
                    
parser.add_argument("-r", "--round",
                    help="Number of digits to round",
                    type=int,
                    required=False)
                    
args = parser.parse_args()

with open(args.file, 'r') as sequence:
    dna = sequence.read().upper()

def validate(dna):
    invalid = re.finditer("[^ATCG]+", dna)
    match = len([invalid])
    if match:
        for error in invalid: 
            print(f"Existen bases invalidas en el archivo introducido: ")
        return(0)
    else:
        return(1)
        
    
