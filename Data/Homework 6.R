# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# putting the install commands as comments allows you to quickly run all of the load
# library commands without accidentally re-installing them all
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



#To get Uniprot t"biogrowth"#To get Uniprot to load correctly 
# BiocManager::install("GenomicAlignments")

Biostrings:writeXStringSet(Bacillus_liceniformis_DNA_sequence.fasta_1.fasta)
a

#Reading the file into R
accession_numbers <- read.csv("UniProt_Entry_Numbers.txt")
accession_numbers$accession.numbers # access the column containing your data
# replace a random space that got into the second value (could edit in the original file, too)
accession_numbers$accession.numbers[2] <- gsub("P00780 ", "P00780", accession_numbers$accession.numbers[2])
#Formating to the list of accession numbers 
# accession_numbers<- c("P29142","P00780","P35835","P04189","P00783")
accession_numbers <- accession_numbers[,1] # another way to access the first column of your data

#Getting the Gene Ontology (GO) terms PROTEIN

AccesionNumbers<-GetProteinGOInfo(accession_numbers)

#Plot results using plotGoInfo()
# output to a folder on your computer, called "PlotGoInfo"
PlotGoInfo(AccesionNumbers, directorypath = "PlotGoInfo")

#find information on pathologies associated with my gene
# save to a variable
GPB <- GetPathology_Biotech(accession_numbers)

#find information on any diseases
# save to variable
diseases <- Get.diseases(GPB)
#Error in UseMethod("select") : 
#no applicable method for 'select' applied to an object of class "character"
# no known associated diseases

#access structural information using the protti package
# save to variable
uniprots <- fetch_uniprot(accession_numbers)

#3D structures for my gene 
# save to variable
af_pred <- fetch_alphafold_prediction(accession_numbers)
