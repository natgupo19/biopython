from os import name
from Bio import Entrez
from pprint import pprint 
import argparse
import re

parser = argparse.ArgumentParser(
    description="")

parser.add_argument("-t", "--termino",
                    help="Introduzca el organismo y los genes en el siguiente formato: 'Organismo1: Gene1, Gene2; Organismo2: Gene3, Gene4...'",
                    type=str,
                    required=True)

parser.add_argument("-o", "--output",
                    help="Archivo de salida",
                    required=False)

args = parser.parse_args()

Entrez.email = "natgupo@lcg.unam.mx"

def create_term(info):
    
    terms = []
    
    for inf in info.split('; '):
        inf = inf.split(':')
        term = inf[0] + "[Orgn] AND ("
        i = True
        for gen in inf[1].replace(' ', '').split(','):
            if i:
                term += gen + "[Gene]"
            else:
                term += " OR " + gen + "[Gene]"
            i = False
        term += ")"
        
        terms += [term] 
        
        
    return(terms)

def term_ids(terms):
    db_ids = []
    for termino in terms:
        organism = termino.split("[Orgn]")
        organism = organism[0]
        handle = Entrez.egquery(term = termino)
        record = Entrez.read(handle)
        print(f"\n{organism}\n")
        
        dbases = []
        for row in record["eGQueryResult"]:
            if row["Count"] != "0" and row["Count"] != "Error":
                dbases += [row["DbName"]]
            
        
        for dba in dbases:
            ids = []
            db_record = Entrez.read(Entrez.esearch(db = dba, term = termino))
            ids = db_record["IdList"]
            dba = dba.split(",")
            db_id = dba[0] + ": "
            db_ids = db_id + str(ids).replace('[', '').replace(']', '').replace("'", "")
            print(f"{db_ids}\n")
    return()

termino = create_term(args.termino)
data_bases = term_ids(termino) 
                
