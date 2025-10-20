#09/08/25
#Función que estima los REPS de una secuencia genómica dada
library("Biostrings")
library("dplyr")

setwd("C:/Users/gusta/Desktop/CursosR ladys/Pp/2012")
REPS <- readDNAStringSet("RepsPass1.fasta")
DNASeq <- readDNAStringSet("GCF_000315235.1_ASM31523v1_genomic.fna")

  ExpCol <- c("IDSeq","IDRep","Start","End","Ident","SeqRes")
  Exp <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(Exp) <- ExpCol
  for (x in 1:length(REPS)) {
    b <- as.character(REPS[x])
    Umbral <- (1-(2/nchar(b)))
    for (y in names(DNASeq)) {
      a <- as.character(DNASeq[y])
      
      for (z in 1:(nchar(a)-nchar(b))) {
        temp <- substr(a,z,(z+nchar(b)-1))
        Valor <- stringdist::stringsim(temp,b, method = "osa")
        if(Valor>= Umbral) {
          
          Exp <- rbind(Exp, data.frame(IDSeq = y, IDRep = b, Start = z, End = (z+nchar(b)-1), Ident = Valor, SeqRes = temp ))
        }
      }
    }
    print(x)
  }
  
  Exp
  write.csv(Exp, "DocReps2012.csv" ,row.names = FALSE)
Salvar <- REPS(REPS = reps, DNASeq = DNASeq)

Salvar 




