library(vegan)
library(ggplot2)
library(tidyverse)


# load the species tables from a full run
species <- read.csv("data/species_matrix.csv", row.names = 1)

# we want to have samples in the rows and species in the columns
species <- t(species)

samples <- read.csv("data/sample_meta.csv", stringsAsFactors = TRUE)


samples$sr <- specnumber(species)
samples$shannon <- diversity(species, index = "shannon")


# check some relationships

boxplot(sr ~ building, data = samples)
boxplot(shannon ~ building, data = samples)

sp <- data.frame(species)
sp$samp <- rownames(sp)

library(reshape)

dat <- melt(species)

ggplot(dat, aes(x = X1, y = value)) +
    geom_bar(stat = "identity", aes(fill = X2))

sp2 <- data.frame(x = colnames(species), y = colSums(species))
ggplot(sp2, aes(x = factor(x), y = log(y), fill = x)) +
    geom_bar(stat = "identity")


###


mds <- metaMDS(species, distance = "bray")
plot(mds, type = "text", display = "sites")
plot(mds)
ordihull(mds, groups = samples$building, draw = "polygon", col = 2:3)


ano <- anosim(species, samples$building, distance = "bray")
ano

out <- c(
    "Salmo salar", "Sus scrofa", "Bos taurus", "Vulpes vulpes",
    "Canis lupus", "Ovis aries", "Gallus gallus", "Homo sapiens",
    "unknown", "Macropus giganteus"
)
sp3 <- species[, !colnames(species) %in% out]

index <- -which(rowSums(sp3) == 0)
sp4 <- sp3[-which(rowSums(sp3) == 0), ]
s2 <- samples[index, ]

s2$sr <- specnumber(sp4)
s2$shannon <- diversity(sp4, index = "shannon")
boxplot(sr ~ building, data = s2)
boxplot(shannon ~ building, data = s2)
t.test(shannon ~ building, data = s2)

mds <- metaMDS((sp4 > 0), distance = "bray")
plot(mds, type = "text", display = "sites")
plot(mds, type = "text")
ordihull(mds, groups = s2$building, draw = "polygon", col = 1:2)

ano <- anosim(sp4, s2$building, distance = "bray")
ano
