'''
NAME
    GenBank_info.py
    
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
    quality_records
'''