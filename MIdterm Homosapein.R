#Ava Ghannam
#Blast results 
#Gene:Homo SPaiens hbb gene for beta globin
#Acccession Number: LC121775
#Length 696

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
# install.packages("seqinr")
library(seqinr)
library(ape)
# install.packages("phangorn")
library(phangorn)
setwd("Data/")
# install.packages("UniProtR")
# install.packages("protti")
library(UniprotR)
library(protti)
library(r3dmol)
# install.packages("biomaRt")
library("GenomicAlignments")
# install.packages("ape")
library(ape)
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
library(seqinr)
library(ape)
library(phangorn)
library(devtools)
library(protti)
library("GenomicAlignments")
library("UniprotR")
library("r3dmol")



#Reading sequence fasta file 
# sequences.fasta not present in directory!
dna_string<- readDNAStringSet("sequences.fasta")

#inputting msa alignment
msaalignment <- msa(dna_string)
msaalignment

#translating sequence to a proein 
# this first translate command wasn't working. The second one did
protein_sequence <- Biostrings::translate(dna_string)
individual6 <- dna_string[[6]] #excellent
amino_acid_sequences <- Biostrings::translate(individual6)
amino_acid_sequences
AA_seq <- as.character(amino_acid_sequences)

#creating an ouput file
output_file <- "individual6.fasta"

writeXStringSet(AAStringSet(AA_seq), file = output_file,
                format = "fasta", width = 60)

#giving the accesion number to a variable 
accession_number <- ("A0A0J9YWK4")

#getting pathologey/ diseases information
GP<-GetPathology_Biotech(accession_number)
GP
Get.diseases(GP, directorypath = getwd())

#3-d structure for protein
FAlpha<-fetch_alphafold_prediction(accession_number)
View(FAlpha)  
