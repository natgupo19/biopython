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
    Programa que identidica las regiones ricas en AT de 
    un archivo con una secuencia de DNA

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
                                        
args = parser.parse_args()

with open(args.file, 'r') as sequence:
    dna = sequence.read().upper()

def validate(dna):
    invalid = re.finditer("[^ATCG+]", dna)
    match = len([*invalid]) 
    if not match:
        return(1)
    if match:
        invalid = re.finditer("[^ATCG+]", dna)
        for error in invalid:
            print (f"\nSe encontaro un caracter invalido: {error.group()} en la posicion {error.span()}\n")
        return(0)
    
def at_regions(dna, at = 13):
    at_rich = re.finditer("([AT]+)", dna)
    match_2 = len([*at_rich])
    if match_2:
        at_rich = re.finditer("([AT]+)", dna)
        for regions in at_rich:
            if len(regions.group()) >= at:
                print(f"\nSe encontro la region rica en AT: {regions.group()} en la posicion {regions.span()}\n")
    else:
        print("\nNo se encontraron regiones ricas en AT\n")
            
if validate(dna):
    if args.size:
        result = at_regions(dna, args.size)
    
    else:
        result = at_regions(dna)
        

    
        
    
            
