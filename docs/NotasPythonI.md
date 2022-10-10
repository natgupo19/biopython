# DICCIONARIOS

Class: BIOPYTHON
Created: August 30, 2022 12:37 PM
Nombre: Natalia Gutiérrez
Type: NOTAS

Ejercicio de Rosalind

```bash
s = "We tried list and we tried dicts also we tried Zen"a = s.split(' ')dicc={}for word in a:  dicc[word]= a.count(word)  for key, value in dicc.items():    #print (key)     #print (value)    print("{} {}".format(key, value))

```

cadena=[]
cadena= input("Introduce el texto que desees, no mas de 10000 letras \n")
lista=cadena.split(' ')
diccionario={}
sinEspacios= cadena.replace(" ","")
if len(sinEspacios) <= 10000:
for word in lista:
cuenta= lista.count(word)
diccionario[word]= cuenta
for word, valor in diccionario.items():
print(word, end=' '),
print(valor)
else:
print("El texto no debe contener más de 10000 letras")

```bash
# METODO KEYS
# lista con el nombre de llaves
Ferecuencia.keys()

# METODO ITEMS
# lista con nombre de llaves y frecuencia de cada una
Frecuencia.items()

# Funcionamiento de for en el diccionario
for key, item in Frecuencia.items()
	print(key)
	print (item)
miau
4
guau
7

# POPITEM
# borrar el último valor 
Cereal = {
	"Marca": "Kellogs

# FROMKEYS
# Devuelve un diccionario con las claves y el valor especificados.

# CLEAR
# Elimina todos los elementos del diccionario

# VALUE
# Devuelve una lista de todos los valores en el diccionario

# SETDEFAULT
# Devuelve el valor de la clave especificada. Si la clave no existe: inserte la clave, con el valor especificado
```