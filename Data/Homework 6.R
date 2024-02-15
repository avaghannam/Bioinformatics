if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)
install.packages("phangorn")
library(phangorn)
setwd("Data/")
install.packages("UniProtR")
install.packages("protti")
library(UniprotR)
library(protti)
library(r3dmol)
install.packages("biomaRt")
library("GenomicAlignments")



#To get Uniprot t"biogrowth"#To get Uniprot to load correctly 
BiocManager::install("GenomicAlignments")

Biostrings:writeXStringSet(Bacillus_liceniformis_DNA_sequence.fasta_1.fasta)
a

#Reading the file into R
read.csv("UniProt_Entry_Numbers.txt")

#Formating to the list of accession numbers 
accession_numbers<- c("P29142","P00780","P35835","P04189","P00783")

#Getting the Gene Ontology (GO) terms PROTEIN
AccesionNumbers<-GetProteinGOInfo(accession_numbers)

#Plot results using plotGoInfo()
PlotGoInfo(AccesionNumbers)

#find information on pathologies associated with my gene
GetPathology_Biotech(accession_numbers)

#find information on any diseases
Get.diseases(accession_numbers)
#Error in UseMethod("select") : 
#no applicable method for 'select' applied to an object of class "character"

#access structural information using the protti package
fetch_uniprot(accession_numbers)

#3D structures for my gene 
fetch_alphafold_prediction(accession_numbers)
