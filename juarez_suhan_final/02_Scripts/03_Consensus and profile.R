###################################
###### Consensus and Profile
####################################

#Ejercicio a resolver:

#A partir de 10 secuencias de 1000 pb en una matriz, generar secuencia consenso y,
#Una cadena de consenso y una matriz de perfiles para la colección.


#Primero cargamos las secuencia de 1000 pb (fueron tomadas del NCBI y cortadas a los primero 1000pb), yo le spuse su nombre común
library(Biostrings)
Bear<-readDNAStringSet("01_Raw_Data/Bear.fasta")
names(Bear)<-c("Bear")
Canis<-readDNAStringSet("01_Raw_Data/Canis.fasta")
names(Canis)<-c("Canis")
Castor<-readDNAStringSet("01_Raw_Data/Castor.fasta")
names(Castor)<-c("Castor")
Zebra<-readDNAStringSet("01_Raw_Data/FishZebra.fasta")
names(Zebra)<-c("ZebraF")
Gallus<-readDNAStringSet("01_Raw_Data/Gallus.fasta")
names(Gallus)<-c("Gallus")
Human<-readDNAStringSet("01_Raw_Data/Human.fasta")
names(Human)<-c("Human")
Monkey<-readDNAStringSet("01_Raw_Data/Monkey.fasta")
names(Monkey)<-c("Monkey")
Mouse<-readDNAStringSet("01_Raw_Data/Mouse.fasta")
names(Mouse)<-c("Mouse")
Perca<-readDNAStringSet("01_Raw_Data/Perca.fasta")
names(Perca)<-c("Perca")
Pig<-readDNAStringSet("01_Raw_Data/Pig.fasta")
names(Pig)<-c("Pig")

#Transformamos cada secuencia en una matriz
Bear_M<-as.matrix(Bear)
Canis_M<-as.matrix(Canis)
Castor_M<-as.matrix(Castor)
Zebra_M<-as.matrix(Zebra)
Gallus_M<-as.matrix(Gallus)
Human_M<-as.matrix(Human)
Monkey_M<-as.matrix(Monkey)
Mouse_M<-as.matrix(Mouse)
Perca_M<-as.matrix(Perca)
Pig_M<-as.matrix(Pig)

##Con la función rbind, combinamos todas las matrices en una sola
TLR1_matrix<-rbind(Bear_M,Canis_M,Castor_M,Zebra_M,Gallus_M,Human_M,Monkey_M,Mouse_M,Perca_M,Pig_M)

##Con esta función verificamos las dimensiones de nuestra matriz
dim(TLR1_matrix)  ##Nos indica una matriz de 10 filas por 1000 columnas

##Vamos a preguntar el número de bases que se repitan en cada columna
TLR1_matrixD<-as.data.frame(TLR1_matrix)
names(TLR1_matrixD)

##vamos a usar el siguiente ciclo while para calcular la matriz de frecuencias de nucleótidos

A_<-c()
T_<-c()
G_<-c()
C_<-c()
i<-1

while (i <= 1000) {
  Columna<-TLR1_matrixD[,i]
  Columna
  
  if (sum((Columna == "A")) > 0) {
    A_<-c(A_,sum(Columna == "A"))
  } else {
    A_<-c(A_,0)
  }
  
  if (sum((Columna == "G")) > 0) {
    G_<-c(G_,sum(Columna == "G"))
  } else {
    G_<-c(G_,0)
  }

  if (sum((Columna == "C")) > 0) {
    C_<-c(C_,sum(Columna == "C"))
  } else {
    C_<-c(C_,0)
  }
  
  if (sum((Columna == "T")) > 0) {
    T_<-c(T_,sum(Columna == "T"))
  } else {
    T_<-c(T_,0)
  }
  
  i<-i+1
  
}

##Colocamos una sola matriz consenso
Frequency_matrix<-rbind(A_,C_,G_,T_)
Frequency_matrix

##Ahora debemos colocar la secuencia consenso con lo siguiente:
Consensus_Seq<-c()
Variantes<-c() ##Si tenemos una base con la misma frecuencia
Base_pos<-c() ##Que nos de la posición de esa variante
i<-1

while (i <= 1000) {
  Col<-Frequency_matrix[,i]
  Reng1<-Frequency_matrix[1,i]
  Reng2<-Frequency_matrix[2,i]
  Reng3<-Frequency_matrix[3,i]
  Reng4<-Frequency_matrix[4,i]
  
  if (max(Col) == Reng1 & max(Col) == Reng2 & max(Col) == Reng3 & max(Col) == Reng4) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"W")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1 & max(Col) == Reng2 & max(Col) == Reng3) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"X")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1 & max(Col) == Reng4 & max(Col) == Reng3) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"Z")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng2 & max(Col) == Reng3 & max(Col) == Reng4) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"Z")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1 & max(Col) == Reng2 & max(Col) == Reng4) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"Y")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1 & max(Col) == Reng2) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"C")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1 & max(Col) == Reng3) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"G")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1 & max(Col) == Reng4) {
    Consensus_Seq<-c(Consensus_Seq,"A")
    Variantes<-c(Variantes,"T")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng2 & max(Col) == Reng3) {
    Consensus_Seq<-c(Consensus_Seq,"C")
    Variantes<-c(Variantes,"G")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng2 & max(Col) == Reng4) {
    Consensus_Seq<-c(Consensus_Seq,"C")
    Variantes<-c(Variantes,"T")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng3 & max(Col) == Reng4) {
    Consensus_Seq<-c(Consensus_Seq,"G")
    Variantes<-c(Variantes,"T")
    Base_pos<-c(Base_pos,i)
  } else if (max(Col) == Reng1) {
    Consensus_Seq<-c(Consensus_Seq,"A")
  } else if (max(Col) == Reng2) {
    Consensus_Seq<-c(Consensus_Seq,"C")
  } else if (max(Col) == Reng3) {
    Consensus_Seq<-c(Consensus_Seq,"G")
  } else {
    Consensus_Seq<-c(Consensus_Seq,"T")
  }  
    
  i<-i+1
  
}

Consensus_Seq
Variantes
Base_pos

###Podemos guardarlo para colocar en su momento las posiciones de variación
names(Base_pos)<-Variantes
Base_pos

##Letras extra del alfabeto de nucleotidos
## X == "C" or "G"
## Z == "G" or "T"
## X == "C" or "T"
## W == "C" or "G" or "T"

###Guardamos la información
##Abrimos la paquetería correspondiente para transformar y guardar nuestros datos de la secuencia de proteínas
##En un formato necesario, como es FASTA
library(gridExtra)
library(datasets)
library(seqinr)
setwd("03_Output/Processed data/") ##Nos posicionamos en la carpeta correspondiente, de salida de datos
write.csv(Frequency_matrix, file = "Matriz_frecuencias.csv")
write.fasta(sequences = Consensus_Seq,names = "Secuencia consenso",file.out="Sequence_consensus.fasta")
write.csv(Base_pos, file = "Matriz_secundaria.csv")
write.fasta(sequences = Base_pos,names = "Secuencia consenso secundaria",file.out="Sustitutciones_consenso.fasta")
dev.off()


###Construyamos un árbol a partir de las secuencias ya empleadas y la consenso, ya tenemos un archivo de esto en datos crudos

TLRs_animals<-readDNAStringSet("01_Raw_Data/TLR1Animals.fasta")
TLRs_animals

Idnomb<-c("Bear","Canis","Castor","Consenso","ZebraFish","Gallus","Homo","Macaca","Mus","Perca","Sus")
names(TLRs_animals)<-Idnomb
TLRs_animals

library(msa)
TLR1<- msa(TLRs_animals, "ClustalW")
print(TLR1, show="alignment")
print(TLR1, show="standardParams")

###Construcción del árbol
TLR1tree<-msaConvert(TLR1, type=c("seqinr::alignment"))

library(seqinr)
dClu<-dist.alignment(TLR1tree, matrix= c("identity"),gap = FALSE) ##genrar matriz de distancias
dClu

library(ape)
TLRsarbol<-nj(dClu)
plot(TLRsarbol, main="Árbol filogenético del receptor TLR1 en distintas especies")

