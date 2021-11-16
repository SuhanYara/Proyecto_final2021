###################################
###### Locating Restriction Sites
####################################

#Ejercicio a resolver:

#A partir de una secuencia de al menos 1000 pb ubicar la posición y longitud de cada palindromo inverso en la cadena
##Palindromo inverso: serán secuencias que tienen la misma secuencia en su complemento(de dercha a izquierda)

#Primero cargamos una secuencia de 1000 pb
library(Biostrings)
Bacillus<-readDNAStringSet("01_Raw_Data/Bacillus_pseudomycoides.fasta")
Bacillus

Bac_matrix<-as.matrix(Bacillus)
Bac_matrix
Sec_Bacillus<-as.vector(Bac_matrix)
Sec_Bacillus
length(Sec_Bacillus)


##Primero necesitamos la cadena complementaria (C-G y A-T), yo usé lo siguiente:
#Busqué la posición de las bases y después un reemplazo por el complementario
Adenin<-which(Sec_Bacillus =="A")
Adenin
Timin<-which(Sec_Bacillus =="T")
Timin
Citocin<-which(Sec_Bacillus =="C")
Citocin
Guanin<-which(Sec_Bacillus =="G")
Guanin

Complement<-replace(Sec_Bacillus,Adenin,"T")
Complement<-replace(Complement,Timin,"A")
Complement<-replace(Complement,Citocin,"G")
Complement<-replace(Complement,Guanin,"C")

#Guardé en un objeto esto:
Complement

#Para hacer una comparación con la cadena principal, reacomodar la cadena complementaria de final a inicio
ReversoC<-Complement[1000:1]
ReversoC

#Ahora hay que buscar los palindromos
Posicion<-1
Palindromo4<-c()

for (Posicion in 1:length(Sec_Bacillus)) {
  if (all(Sec_Bacillus[Posicion:(Posicion+3)]==ReversoC[Posicion:(Posicion+3)])) {
    Palindromo4<-c(Palindromo4, Posicion)
  }
  Posicion<-Posicion +1
}

##Aquí vemos el número de palindromos
Palindromo4_Número<-length(Palindromo4)
Palindromo4_Número

##Aquí las posiciones
Posiciónes_Palindromos4<-Palindromo4
Posiciónes_Palindromos4

##Con esto veo la secuencia de cada palindromo
Sec_Bacillus[117:120]
Sec_Bacillus[881:884]


###Si sólo existen un límitado número de palindromos de tamaño 4, como el ejemplo, sólo podríamos buscar
#en esas posiciones algunos con tamaño mayor. O aplicar el comando anterior y adicionar para 5 y hasta 12

##Aquí únicamente haré una comparación simple, porque sólo tengo 2

Posiciónes_Palindromos4

##Marca 117 y 881, y esas buscaré con una comparación simple para ver si son iguales

##Para palindromos de 5
all(Sec_Bacillus[117:121]==ReversoC[117:121])
all(Sec_Bacillus[881:885]==ReversoC[881:885])

##Ambas son falsas, así que no podremos continuar si sólo tenemos estos 







