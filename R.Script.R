ds=5
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
#Set Working directory for Data and Homework folder ####

getwd()
#Checked that directory was in place ####
mySequences01 <- readDNAStringSet("Bacillus_licheniformis_DNA_sequence.fasta_1.fasta")
mySequences02 <- readDNAStringSet("Bacillus_licheniformis_DNA_sequence.fasta_2.fasta")
mySequences03 <- readDNAStringSet("Bacillus_licheniformis_DNA_sequence.fasta_3.fasta")
mySequences04 <- readDNAStringSet("Bacillus_licheniformis_DNA_sequence.fasta_4.fasta")
mySequences05 <- readDNAStringSet("Bacillus_licheniformis_DNA_sequence.fasta_5.fasta")
mycombinedSeq <- readDNAStringSet("msa.txt")

mySequenceFile <- c(mySequences01, mySequences02, mySequences03, mySequences04, mySequences05)

mySequenceFile 
#### Creating a consensus ####
alignment_set <- DNAStringSet(mySequenceFile)
consensus <- consensusString(alignment_set)
print(consensus)

#### Creating an msa alignment ####
myFirstAlignment <- msa(mySequenceFile)
myFirstAlignment

print(myFirstAlignment, show="complete")

####Function to count gaps in a consensus ####
count_gaps <- function(sequence) 
  #Search for the presence of "-" in the sequence
  gap_positions <- grepl("-", sequence)

#Count the number of TRUE values (i.e., gaps)
num_gaps <- sum("gap_positions")

return(num_gaps)

#Count gaps in each sequence of the consensus
num_gaps_in_sequences <- sapply(consensus, count_gaps)

#Total gaps in the alignment
total_gaps <- sum(num_gaps_in_sequences)


#Print the total number of gaps
print(total_gaps)
#1 gap

#####Calculate the width of the alignments ####
alignment_length <- width(consensus)

#Print the length of the alignments
print(alignment_length)
#ALignmnet length is 1950

####Calculate GC content ####
# Convert the alignment to a DNAStringSet object
alignment <- DNAStringSet(consensus)

# Calculate the GC content for each position in the alignment
c_content <- vcountPattern("C", alignment) / width(alignment)
G_content <- vcountPattern("G", alignment) / width(alignment)

# Print the GC content for each position
print(c_content)
print(G_content)
#[1] 0.2758974
#[1] 0.2071795
#GC content is .482 = 48.2%

####Convert alignment to SeqinR format ####
Bacilus_liceniformisCom<- msa(mycombinedSeq)
Bacilus_liceniformisCom

Bacilus_liceniformisCom <- msaConvert(Bacilus_liceniformisCom, type="seqinr::alignment")
# Compute distance matrix using SeqinR functions
# Assuming you want to compute identity distance
d <- dist.alignment(Bacilus_liceniformisCom)

# View the distance matrix
print(d)

####Creat Phylogenetic Tree ####
Env2Tree <- nj(d)
plot(Env2Tree, main="Phylogenetic Tree of Bacilus_liceniformis Gene Sequences")

#### Translate DNA sequence into AA sequence ####
dna_sequences <- readDNAStringSet("Bacillus_licheniformis_DNA_sequence.fasta_1.fasta")
amino_acid_sequences <- translate(dna_sequences)amino_acid_sequences

####Convert alignment to phangorn ####
Alignment_phyDat <- msaConvert(myFirstAlignment, type="phangorn::phyDat")
Alignment_phyDat
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")



a
a
