source("dirForDEanalysis.R")
library(tximport, lib.loc = rLibPackageDir)
library(rhdf5, lib.loc = rLibPackageDir)
library(DESeq2, lib.loc = rLibPackageDir)



#read in full sample names:
sampleNameFull <- read.csv(sampleNameFullFile, header = FALSE)

#read in phenotype file:
pheno <- read.csv(phenoDir, header = TRUE, sep="\t")

#create sample files array:
##remember, I lauched R in the same directory as this file:
files <- file.path(sampleDir, sampleNameFull[,1], "abundance.tsv")
names(files) <- pheno$sampleName

#importing the abundance.tsv data:
txi.tx <- tximport(files, type = "kallisto", txOut = TRUE)

#creating the DESeq object for analysis
dds <- DESeqDataSetFromTximport(txi.tx, pheno, ~sex)

#filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#DE analysis
dds <- DESeq(dds)
res <- results(dds)

#reordering matrix by pvalue from smallest to largest
resOrdered <- res[order(res$pvalue),]

#writing the results to file
date <-format(Sys.time(), format = "%Y-%m-%d_%H_%M_%S")
outputDir <- file.path(DESeqOutputDir, date)
dir.create(outputDir)
outputFile <- file.path(outputDir, "DESeq2-MaleVsFemale.txt")
write.csv(resOrdered, outputFile)

