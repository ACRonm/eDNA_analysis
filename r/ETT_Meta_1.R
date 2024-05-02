# if you want to run this on your own computer you need to uncomment the code below to install all packages
# be prepared that may take a while...

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("dada2")
# BiocManager::install("phyloseq")
# BiocManager::install("microbiome")
# BiocManager::install("DECIPHER")
# BiocManager::install("msa")

# install.packages("tidyverse")
# install.packages("phangorn")


# Load all required packages
sapply(c("dada2", "phyloseq", "tidyverse", "microbiome", "msa", "DECIPHER", "phangorn", "ape", "seqinr"),
    require,
    character.only = TRUE
)


### load and examine single fasta sequences###
sequencefile <- "../data/12S_refdb_curated.fasta"
mySequences <- readAAStringSet(sequencefile)
mySequences

### align the sequences#
# now you can align your individual sequences using the msa package
alignment <- msa(mySequences)
print(alignment, show = "complete")

### there are many things one can do with the alignment. for now we will convert our msa alignment for use in a different package###
biosec <- msaConvert(alignment, type = "seqinr::alignment")


### Now, we will create a distance matrix to determine how related these sequences are from each other by using dist.alignment from the seqinr package. These functions compute a matrix of pairwise distances from aligned sequences using similarity (Fitch matrix, for protein sequences only) or identity matrix (for protein and DNA sequences).The resulting matrix contains the squared root of the pairwise distances. ###

d <- dist.alignment(biosec, "identity")

heatmap(as.matrix(d))
# maybe need to change lable sizes...
heatmap(as.matrix(d), cexRow = 0.5, cexCol = 0.5)
# save as png

png("ETT_single_sequences_heatmap.png", width = 10, height = 10, units = "in", res = 600)
heatmap(as.matrix(d), cexRow = 0.5, cexCol = 0.5)
dev.off() # Close the png device

### the heatmap shows you how many similarity groups are in your alignment. values close to 0 show that sequences have no differences between one another. values close to 1 sho highly disssimilar sequences

### now, create a neighbor joined tree from your distance matrix###

tree <- nj(d)

# save the tree as a pdf
pdf("ETT_single_sequences_tree.pdf")
plot(tree, type = "fan", cex = 0.5)
dev.off()

png("ETT_single_sequences_tree.png", width = 10, height = 10, units = "in", res = 300)
plot(tree, type = "fan", cex = 0.5)
dev.off()
