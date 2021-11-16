####################################
###### 1. Transcribing DNA into RNA 
####################################

#Ejercicio a resolver:

#A partir de una secuencia de al menos 1000 pb
#Transcribir de ADN a ARN

#El problema implica sustituir la Timina, sólo en el ADN por Uracilo en ARN.

#Vamos a empezar con una cadena pequeña que contenga las cuatro bases
library(Biostrings)
A<-DNAString("AGCT")
A
#Ahora la transformamos en un vector para trabajar con ella, ya que sino se leerá como cadena de ADN del paquete 
A<-as.vector(A)
class(A)  #Revisamos sus características
length(A) #Tamaño
A

#Aquí hay una característica importante, ya que la secuencia se tomó desde un programa base
#ya tiene niveles definidos...
levels(A)

##si deseamos reemplzar T por U, éste último no se encuentra en los niveles y no lo reemplzará
##Por ello vamos a agregar ese nivel con el siguiente código
levels(A)<-c(levels(A),"U")
levels(A)

##Ya que lo tenemos agregado, vamos a usar which para conocer la posición de las T en el vector
Tpos<-which(A=="T")
Tpos
##Y lo guardamos en un objeto

##Ahora simplmente usamos la función replace para sutituirlo y listo
replace(A,Tpos,"U")


#Ahora vamos a usar una secuencia de nucleotidos bajada del NCBI, de una especie de Bacillus
#con los 1000 nucleotidos
#En secuencias grandes, por el tipo de datos, no se puede aplicar lo anterior, 
#por ello, se puede unar una simple función de sustitución
SeqBacillus<-readDNAStringSet("01_Raw_Data/Bacillus_pseudomycoides.fasta")
SeqBacillus

#Aquí se encuentra el siguiente código
RNABacillus<-gsub("T","U",SeqBacillus)
RNABacillus

####################################
###### 2. Transcribing RNA into protein
####################################

#Ejercicio a resolver:

#A partir de una secuencia de ARN Transcribir a proteína
#Usamos la siguiente secuencia

Mirounga<-readDNAStringSet("01_Raw_Data/Mirounga_MYD88.fasta")
Mirounga<-as.matrix(Mirounga)
Mirounga<-as.vector(Mirounga)
Mirounga


##Según lo anterior, lo transformamos a RNA
RNAMirounga<-gsub("T","U",Mirounga)
RNAMirounga

###NOTRNAMirounga###NOTA: Se transform los datos para poder acceder y trabajar con ellos

##Haremos VECTORES que corresponden a los codones para poder trabajar con los datos de la secuencia y que se comparen
F1<-c("U","U","U")
F2<-c("U","U","C")
L1<-c("U","U","A")
L2<-c("U","U","G")
L3<-c("C","U","U")
L4<-c("C","U","C")
L5<-c("C","U","A")
L6<-c("C","U","G")
I1<-c("A","U","U")
I2<-c("A","U","C")
I3<-c("A","U","A")
M1<-c("A","U","G")
V1<-c("G","U","U")
V2<-c("C","C","C")
V3<-c("G","U","A")
V4<-c("G","U","G")
S1<-c("U","C","U")
S2<-c("G","U","C")
S3<-c("U","C","A")
S4<-c("U","C","G")
P1<-c("C","C","U")
P2<-c("C","C","C")
P3<-c("C","C","A")
P4<-c("C","C","G")
T1<-c("A","C","U")
T2<-c("A","C","C")
T3<-c("A","C","A")
T4<-c("A","C","G")
A1<-c("G","C","U")
A2<-c("G","C","C")
A3<-c("G","C","A")
A4<-c("G","C","G")
Y1<-c("U","A","U")
Y2<-c("U","A","C")
STOP1<-c("U","A","A")
STOP2<-c("U","A","G")
STOP3<-c("U","G","A")
H1<-c("C","A","U")
H2<-c("C","A","C")
Q1<-c("C","A","A")
Q2<-c("C","A","G")
N1<-c("A","A","U")
N2<-c("A","A","C")
K1<-c("A","A","A")
K2<-c("A","A","G")
D1<-c("G","A","U")
D2<-c("G","A","C")
E1<-c("G","A","A")
E2<-c("G","A","G")
C1<-c("U","G","U")
C2<-c("U","A","C")
W1<-c("U","G","G")
R1<-c("C","G","U")
R2<-c("C","G","C")
R3<-c("C","G","A")
R4<-c("C","G","G")
R5<-c("A","G","A")
R6<-c("A","G","G")
S1<-c("A","G","U")
S2<-c("A","G","C")
G1<-c("G","G","U")
G2<-c("G","G","C")
G3<-c("G","G","A")
G4<-c("G","G","G")

###
#Fragmentamos RNAMirounga para que cada codón sea leído y reemplazado por el aa correspondiente

Frag_RNAMirounga<-split(RNAMirounga, ceiling(seq_along(RNAMirounga)/3))

##Creamos un vector vacío para que se guarde la sec de proteínas
ProtMirounga<-c()

##Empleamos el siguiente Script para determinar la secuencia de aminácidos:

for (codon in Frag_RNAMirounga) {
  
  if (all(codon == F1) | all(codon == F2)) {
    ProtMirounga<-c(ProtMirounga, "F")
  } else if (all(codon == L1) | all(codon == L2) | all(codon == L3) | all(codon == L4) | all(codon == L5) | all(codon == L6)) {
    ProtMirounga<-c(ProtMirounga, "L")
  } else if (all(codon == I1) | all(codon == I2) | all(codon == I3)) {
    ProtMirounga<-c(ProtMirounga, "I") 
  } else if (all(codon == M1)) {
    ProtMirounga<-c(ProtMirounga, "M")
  } else if (all(codon == V1) | all(codon == V2) | all(codon == V3) | all(codon == V4)) {
    ProtMirounga<-c(ProtMirounga, "V")
  } else if (all(codon == S1) | all(codon == S2) | all(codon == S3) | all(codon == S4)) {
    ProtMirounga<-c(ProtMirounga, "S")
  } else if (all(codon == P1) | all(codon == P2) | all(codon == P3) | all(codon == P4)) {
    ProtMirounga<-c(ProtMirounga, "P")
  } else if (all(codon == T1) | all(codon == T2) | all(codon == T3) | all(codon == T4)) {
    ProtMirounga<-c(ProtMirounga, "T")
  } else if (all(codon == A1) | all(codon == A2) | all(codon == A3) | all(codon == A4)) {
    ProtMirounga<-c(ProtMirounga, "A")
  } else if (all(codon == Y1) | all(codon == Y2)) {
    ProtMirounga<-c(ProtMirounga, "Y")
  } else if (all(codon == STOP1) | all(codon == STOP2) | all(codon == STOP3)) {
    ProtMirounga<-c(ProtMirounga, "STOP")
  } else if (all(codon == H1) | all(codon == H2)) {
    ProtMirounga<-c(ProtMirounga, "H")
  } else if (all(codon == Q1) | all(codon == Q2)) {
    ProtMirounga<-c(ProtMirounga, "Q")
  } else if (all(codon == N1) | all(codon == N2)) {
    ProtMirounga<-c(ProtMirounga, "N")
  } else if (all(codon == K1) | all(codon == K2)) {
    ProtMirounga<-c(ProtMirounga, "K")
  } else if (all(codon == D1) | all(codon == D2)) {
    ProtMirounga<-c(ProtMirounga, "D")
  } else if (all(codon == E1) | all(codon == E2)) {
    ProtMirounga<-c(ProtMirounga, "E")
  } else if (all(codon == C1) | all(codon == C2)) {
    ProtMirounga<-c(ProtMirounga, "C")
  } else if (all(codon == W1)) {
    ProtMirounga<-c(ProtMirounga, "W")
  } else if (all(codon == R1) | all(codon == R2) | all(codon == R3) | all(codon == R4) | all(codon == R5) | all(codon == R6)) {
    ProtMirounga<-c(ProtMirounga, "R")
  } else if (all(codon == S1) | all(codon == S2)) {
    ProtMirounga<-c(ProtMirounga, "S")
  } else if (all(codon == G1) | all(codon == G2) | all(codon == G3) | all(codon == G4)) {
    ProtMirounga<-c(ProtMirounga, "G")
  } else {
    ProtMirounga<-c(ProtMirounga, codon)
  }
  
}

##Aquí podemos visualizar la secuencia de bases
ProtMirounga

##Abrimos la paquetería correspondiente para transformar y guardar nuestros datos de la secuencia de proteínas
##En un formato necesario, como es FASTA

library(gridExtra)
library(datasets)
library(seqinr)
setwd("03_Output/Processed data/") ##Nos posicionamos en la carpeta correspondiente, de salida de datos
write.fasta(sequences = ProtMirounga,names = "Mirounga_MyD88",file.out="MyD88_aa_NES.fasta")
dev.off()

