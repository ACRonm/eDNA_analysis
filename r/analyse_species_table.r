library(vegan)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggplot2)


# load the species tables from a full run
species <- read.csv("/data/species_matrix_transposed.csv", row.names = 1)

# we want to have samples in the rows and species in the columns
# species <- t(species)  # Comment out this line

samples <- read.csv("data/sample_meta.csv", stringsAsFactors = TRUE)

print(species)

print(samples)

samples$sr <- specnumber(species)
samples$shannon <- diversity(species, index = "shannon")

# plot the abundance of species in each sample from species table
sp2 <- data.frame(x = colnames(species), y = colSums(species))
ggplot(sp2, aes(x = factor(x), y = log(y), fill = x)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Species") +
    ylab("Total Abundance (log(y))")

# save the plot
ggsave("data/total_abundance.png", plot = last_plot(), width = 10, height = 10, dpi = 300)




# check some relationships

boxplot(sr ~ site, data = samples)

boxplot(shannon ~ site, data = samples)

# change x tick angle
ggplot(samples, aes(x = site, y = sr)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_blank()) +
    xlab("Molonglo River")

# save the plot

ggsave("data/sr_boxplot.png", plot = last_plot(), width = 10, height = 10, dpi = 300)
# change x tick angle

ggplot(samples, aes(x = site, y = shannon)) +
    geom_boxplot() +
    theme(axis.text.x = element_text()) +
    xlab("Molonglo River site proximity") +
    ylab("Genetic Diversity (shannon index)") +
    theme(axis.title = element_text(size = 14)) +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.ticks = element_line(size = 1)) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank())

# save the plot

ggsave("data/shannon_boxplot.png", plot = last_plot(), width = 10, height = 10, dpi = 300)


sp <- data.frame(species)
sp$samp <- rownames(sp)

library(reshape)

# Melt the species dataframe
dat <- melt(species)

# Create the bar plot
ggplot(dat, aes(x = variable, y = value)) +
    geom_bar(stat = "identity", aes(fill = value))

sp2 <- data.frame(x = colnames(species), y = colSums(species))
ggplot(sp2, aes(x = factor(x), y = log(y), fill = x)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Species") +
    ylab("Total Abundance (log(y))")

# save the plot
ggsave("data/total_abundance.png", plot = last_plot(), width = 10, height = 10, dpi = 300)


# close plot
dev.off()
###


mds <- metaMDS(species, distance = "bray")

# Create a PNG file
png("data/mds_plot.png", width = 10, height = 10, units = "in", res = 300)

# Plot the NMDS results
plot(mds, type = "text", display = "sites")

# Draw the hulls
ordihull(mds, groups = samples$site, draw = "polygon", col = 2:3)

# Close the PNG device
dev.off()

print(species, samples$site)

print(unique(samples$site))

ano <- anosim(species, samples$site, distance = "bray")
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
