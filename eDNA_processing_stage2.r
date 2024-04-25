# libraries --------------------------------------------------------------------
library(tidyverse)
library(dada2)
library(phyloseq)
library(metacoder)
library(gridExtra)

# set file paths ---------------------------------------------------------------
linux <- FALSE
path_to_fastq <- "./data/ETT_12S_fasta/"
path_filtered_fastq <- "./data/ETT_filtered_sequences" # folder name for new filtered fastq

# fastq files ------------------------------------------------------------------

# forward sequences
fnFs <- sort(list.files(path_to_fastq,
    pattern = "_L001_R1_001.fastq.gz", full.names = TRUE
))
# reverse sequences
fnRs <- sort(list.files(path_to_fastq,
    pattern = "_L001_R2_001.fastq.gz", full.names = TRUE
))

# sample names
sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)

# trim "ETT" from front of sample names
sample.names <- gsub("ETT_", "", sample.names)

# trim _S* from end of sample names
sample.names <- gsub("_S[0-9]*", "", sample.names)


# create files for filtered reads
filtFs <- file.path(
    path_filtered_fastq,
    paste0(sample.names, "_F_filtered.fastq.gz")
)
filtRs <- file.path(
    path_filtered_fastq,
    paste0(sample.names, "_R_filtered.fastq.gz")
)

# plot quality -----------------------------------------------------------------

# first two fasta files
plotQualityProfile(fnFs[1:2]) # forward

# save plot
ggsave("./data/quality_profile_forward.png", plotQualityProfile(fnFs[1:2]),
    units = "cm", width = 21, height = 15, dpi = 600
)

plotQualityProfile(fnRs[1:2]) # reverse
# save plot
ggsave("./data/quality_profile_reverse.png", plotQualityProfile(fnRs[1:2]),
    units = "cm", width = 21, height = 15, dpi = 600
)

# # filter and trim --------------------------------------------------------------
# system.time(out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#     truncLen = c(140, 140),
#     trimLeft = c(18, 18), # based on V3V4 primers
#     maxN = 0, # max Ns
#     maxEE = c(3, 4), # expected errors
#     truncQ = 2, #
#     rm.phix = TRUE,
#     compress = TRUE,
#     multithread = linux
# ))

# windows can't support multi-thread
# this takes approx 3 minutes
# errors -----------------------------------------------------------------------
errF <- learnErrors(filtFs, multithread = linux)
errR <- learnErrors(filtRs, multithread = linux)

# plot errors ------------------------------------------------------------------
plotErrors(errF, nominalQ = TRUE)
ggsave("./data/errF.png", plotErrors(errF, nominalQ = TRUE),
    units = "cm", width = 21, height = 15, dpi = 600
)


plotErrors(errR, nominalQ = TRUE)
ggsave("./data/errR.png", plotErrors(errR, nominalQ = TRUE),
    units = "cm", width = 21, height = 15, dpi = 600
)

# dereplication ----------------------------------------------------------------
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# adding sample names to the dereplicated reads
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# sequence inference - denoised data -------------------------------------------
dadaFs <- dada(derepFs, err = errF, multithread = linux)
dadaRs <- dada(derepRs, err = errR, multithread = linux)

# denoised data
dadaFs[[1]]
dadaRs[[1]]

# merge paired reads -----------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# ASV table --------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
table(nchar(getSequences(seqtab)))


# remove chimeras --------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab,
    method = "consensus",
    multithread = TRUE, verbose = TRUE
)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
# track reads ------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(
    out, sapply(dadaFs, getN), sapply(dadaRs, getN),
    sapply(mergers, getN), rowSums(seqtab.nochim)
)

colnames(track) <- c(
    "input", "filtered", "denoisedF",
    "denoisedR", "merged", "nonchim"
)
rownames(track) <- sample.names

track
track %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    mutate(building = substr(sample_id, 1, 2)) %>%
    group_by(building) %>%
    summarise_if(is.numeric, sum)

# taxa assign ------------------------------------------------------------------
system.time(mytaxa <- assignTaxonomy(
    seqtab.nochim,
    "./refdb/reference_16S_taxa.fa"
))
dim(mytaxa)

unname(mytaxa[, 1:5])
unname(mytaxa[is.na(mytaxa[, 7]), 1:7])
unname(which(is.na(mytaxa[, 7])))


# asv and taxa -----------------------------------------------------------------
seqtab.nochim.df <- t(seqtab.nochim) %>%
    as.data.frame() %>%
    rownames_to_column("sequence") %>%
    mutate(asv_id = paste0("ASV", 1:nrow(.)))

mytaxa.df <- mytaxa %>%
    as.data.frame() %>%
    rownames_to_column("sequence")

fx <- function(x) ifelse(x < 5, 0, x)
asvcounts <- left_join(mytaxa.df, seqtab.nochim.df) %>%
    relocate(asv_id) %>%
    rowwise() %>%
    mutate_if(is.numeric, fx)

# column name meta data --------------------------------------------------------

meta <- data.frame(s = colnames(asvcounts)) %>%
    filter(grepl("C1", s) | grepl("C2", s)) %>%
    mutate(
        building = substr(s, 1, 2),
        session = substr(s, 3, 3),
        sample = substr(s, 4, 4)
    ) %>%
    dplyr::rename(sample_id = s) %>%
    as_tibble()

# metacoder lineage ------------------------------------------------------------
asvspecies <- asvcounts %>%
    select(asv_id:Species) %>%
    mutate(QueryID = asv_id)
lineage <- asvspecies %>%
    # filter(Phylum == 'Chordata') %>%
    mutate(lineage = paste0(
        "r__Root;", "p__",
        Phylum, ";c__", Class, ";o__", Order,
        ";f__", Family, ";g__", Genus
    )) %>%
    select(QueryID, lineage)

# lineage and counts -----------------------------------------------------------
asvLineage <- left_join(asvcounts, lineage, by = c("asv_id" = "QueryID")) %>%
    select(c("asv_id", "lineage", meta$sample_id)) %>%
    relocate(asv_id, lineage) %>%
    tibble() %>%
    filter(complete.cases(lineage))

# metacoder prep - interactive -------------------------------------------------

parse_tax_data(asvLineage,
    class_cols = "lineage", class_sep = ";",
    class_regex = "^(.+)__(.+)$",
    class_key = c(
        tax_rank = "info",
        tax_name = "taxon_name"
    )
)


obj <- parse_tax_data(asvLineage,
    class_cols = "lineage", # the column that contains taxonomic information
    class_sep = ";", # The character used to separate taxa in the classification
    class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
    class_key = c(
        tax_rank = "info", # A key describing each regex capture group
        tax_name = "taxon_name"
    )
)


obj$data$tax_data <- zero_low_counts(obj, dataset = "tax_data", min_count = 5) # change
no_reads <- rowSums(obj$data$tax_data[, meta$sample_id]) == 0
sum(no_reads)


obj <- filter_obs(obj, target = "tax_data", !no_reads, drop_taxa = TRUE)
print(obj)
obj$data$tax_data <- calc_obs_props(obj, "tax_data")

obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
    cols = meta$sample_id
)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund",
    groups = meta$building, # groups for meta coder (can change to session)
    cols = meta$sample_id
)


# plot metacoder ---------------------------------------------------------------

set.seed(1)
tree <- heat_tree(obj,
    node_label = taxon_names,
    node_size = n_obs,
    node_color = "#green"
    node_size_axis_label = "asv count",
    node_color_axis_label = "Samples with reads",
    layout = "davidson-harel", # The primary layout algorithm
    initial_layout = "reingold-tilford"
) # The layout algorithm that initializes node locations
tree

# save as png
ggsave("./data/heat_tree.png", tree,
    units = "cm", width = 21, height = 15, dpi = 600
)


# species table ----------------------------------------------------------------
speciesTbl <- asvcounts %>%
    select(Species:Site_22) %>%
    mutate(Species = ifelse(is.na(Species), "unknown", Species)) %>%
    group_by(Species) %>%
    summarise_if(is.numeric, sum) %>%
    column_to_rownames(var = "Species") %>%
    as.matrix()

speciesTbl2 <- speciesTbl[rowSums(speciesTbl) != 0, ]
rowSums(speciesTbl2)
t(speciesTbl)

write.csv(speciesTbl2, "./data/species_matrix.csv")
write.csv(meta, "./data/sample_meta.csv", row.names = FALSE)
