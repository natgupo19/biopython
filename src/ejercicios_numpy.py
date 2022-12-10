'''
NAME
    Numpy.py
    
VERSION
    1.0
    
AUTHOR
    Natalia Gutierrez Ponce
    
GITHUB
    https://github.com/natgupo19/biopython
    
DESCRIPTION
    

CATEGORY 
    BioPython
    
USAGE
    py .\src\ejercicios_numpy.py
'''

import numpy as np

# EJERCICIO 1 (parte 1)

# Pasar a array el archivo
array_gMl = np.genfromtxt("data/prod_gml.csv", delimiter = ",")

# Convertir de gramos/mililitro a gramos/litro
array_gL = array_gMl * 1000

# Imprimir el resultado
print(f"\nUnidades en g/L:\n{str(array_gL).replace('[',' ').replace(']',' ')}\n")



# EJERCICIO 1 (parte 2)

# Pasar a array el archivo
costos = np.genfromtxt("data/ind_cost.csv", delimiter = ",")

# Convertir los costos
costos_30 = costos * 1.75
costos_35 = costos * 0.8

# Imprimir los resultados
print(f"\nCostos a 30 grados: {str(costos_30).replace('[',' ').replace(']',' ')}")
print(f"Costos a 35 grados: {str(costos_35).replace('[',' ').replace(']',' ')}\n")



# EJERCICIO 2

# Sacar los costos de 1 gramo/litro de cada condicion
costos_30_gL = costos_30 / array_gL[:,0]
costos_35_gL = costos_35 / array_gL[:,1]



# EJERCICIO 3

# Crear un array de numeros que representen los genes
genes = np.array([1, 2, 3, 4])

# Sacar la diferencia de costos por condicion
dif_costos = costos_30_gL - costos_35_gL

# Imprimir los genes que son mas baratos en cada condicion
print(f"\nGenes mas baratos a 30 C: {genes[dif_costos < 0]}")
print(f"Genes mas baratos a 35 C: {genes[dif_costos > 0]}\n")