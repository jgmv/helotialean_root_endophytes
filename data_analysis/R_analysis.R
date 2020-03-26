### load libraries and functions
library(ape)
library(dplyr)
library(pegas)
library(phytools)
library(plyr)
library(reshape)
source("R_code/functions.R")


### create output folder
if (!dir.exists("output")) dir.create("output")


### load phylogenetic trees
tree_mb <- read.tree("data/tree_bayes.tre")
tree_ml <- read.tree("data/tree_raxml.tre")

# drop repeated tips
tree_mb <- drop.tip(tree_mb, c("Cadophora_antarctica_FMR16056",
  "Mycochaetophora_gentianae"))
tree_ml <- drop.tip(tree_ml, c("Cadophora_antarctica_FMR16056",
  "Mycochaetophora_gentianae"))


### load strains' data
strains <- read.csv("data/strains.csv", h = T, sep = ";", row.names = 1)
refs <- read.csv("data/references.csv", h = T, sep = ";", row.names = 1)
refs <- get_country(refs)

# calculate new variables
strains$con_ratio <- strains$con_len / strains$con_wid
strains$con_vol   <- (4 / 3) * pi * strains$con_len * strains$con_wid^2 


### plot phylogenetic trees with source information
plot_fancy_tree(tree_mb, file = "output/Bayes_tree.pdf")
plot_fancy_tree(tree_ml, file = "output/RAxML_tree.pdf")

# plot map for legend
map_legend(file = "output/map_legend.pdf")


### compare phylogeny with OTU clustering
# subset tree to keep only isolates
tree_strains <- drop.tip(tree_mb, grep("^P[0-9]", tree_mb$tip.label,
  invert = T))
tree_strains <- ladderize(tree_strains, right = F)

# compare phylogeny with OTU clusters
phylo_vs_OTUs(tree_strains, file = "output/tree_OTUs.pdf")


### PCA analysis of morphological traits
# calculate PCAs
pcam <- trait_pca(file = "output/PCA_trait.pdf", include_qual = F)
pcam2 <- trait_pca(file = "output/PCA_trait_cat.pdf", include_qual = T)
trait_pca_media(file = "output/PCA_trait_media.pdf", include_qual = F)


# Mantel test: character vs phylogenetic distances
mantel_morpho(tree_mb, pcam)


### Ancetral State Reconstruction (ASR) analysis
# plot trees with quantitative characters
plot_tree_cont(sqrt(strains$ncon), file = "output/acr_ncon.pdf")
plot_tree_cont(strains$con_len, file = "output/acr_con_len.pdf")
plot_tree_cont(strains$con_wid, file = "output/acr_con_wid.pdf")
plot_tree_cont(sqrt(strains$con_vol), file = "output/acr_con_vol.pdf")
plot_tree_cont(strains$con_ratio, file = "output/acr_con_ratio.pdf")
plot_tree_cont(strains$cma_diam, file = "output/acr_cma_diam.pdf")
plot_tree_cont(strains$pda_diam, file = "output/acr_pda_diam.pdf")
plot_tree_cont(strains$mea_diam, file = "output/acr_mea_diam.pdf")
plot_tree_cont(strains$CMA_pigm, file = "output/acr_cma_pigm.pdf")
plot_tree_cont(strains$PDA_pigm, file = "output/acr_pda_pigm.pdf")
plot_tree_cont(strains$MEA_pigm, file = "output/acr_mea_pigm.pdf")

# load discrete characters data
species_ecology <- read.csv("data/species_ecology.csv", h = T, sep = ";")

# subset tree with one tip per species
tree_species <- species_ASR_tree(tree_mb, species_ecology)

# plot trees with ASR analysis
plot_ASR(file ="output/ASR.pdf", file_reduced = "output/ASR_reduced.pdf")


### calculate phylogenetic signal for quantitative characters
calc_phylosignal(tree_mb, file = "output/phylosig_traits.csv")


### estimate species distribution ranges
# input BLAST results
data <- read.csv("data/blast_nt.csv", h = T, sep = ";")
data_sra <- read.csv("data/blast_sra.csv", h = T, sep = ";")

# plot distribution maps per isolate
#for(i in unique(data$qseqid)) plot_distribution(i)

# plot distribution maps per species
for(i in unique(strains$species)) plot_species_distribution(i)

# get species ranges at different percent identities
strains$range_100 <- get_species_ranges(pident = 100)
strains$range_99 <- get_species_ranges(pident = 99)
strains$range_98 <- get_species_ranges(pident = 98)

# plot species ranges
plot_spp_ranges(file = "output/ranges.pdf")


### relationship of species distribution with ecological variables
# get phylogenetic distance
dist_phylo <- cophenetic(tree_strains)
dist_phylo <- dist_phylo[rownames(strains), rownames(strains)]
var_phylo <- cmdscale(dist_phylo, k = 6)

# get morphological distance
var_morph <- pcam$ind$coord
var_morph <- var_morph[rownames(dist_phylo), ]

# get climatic distance
ecol_var <- read.csv("data/ecological_variables.csv", h = T, sep = ";",
  row.names = 1)
ecol_var <- ecol_var[rownames(dist_phylo), ]
var_clim <- ecol_var[, grep("^bio", colnames(ecol_var))]
dist_clim <- dist(scale(var_clim), diag = T, upper = T)
dist_clim <- as.matrix(dist_clim)
dist_clim <- dist_clim[rownames(dist_phylo), rownames(dist_phylo)]

# get soil distance
var_soil <- ecol_var[, c("pH", "N", "organic_mater", "lime", "avP", "P",
  "avK", "S", "cond", "Na", "K", "Mg", "Ca")]
var_soil <- var_soil[rownames(var_soil) != "P1443" &
  rownames(var_soil) != "P1909", ]
dist_soil <- dist(scale(var_soil), diag = T, upper = T)
dist_soil <- as.matrix(dist_soil)
dist_soil <- dist_soil[rownames(dist_phylo), rownames(dist_phylo)]

# get geographical variables
library(sp)
dist_geo <- spDists(as.matrix(strains[, c("longitude", "latitude")]),
  longlat = T)
rownames(dist_geo) <- rownames(strains)
colnames(dist_geo) <- rownames(strains)
var_geo <- pcnm(dist_geo)$vectors
rownames(var_geo) <- rownames(strains)
var_geo <- var_geo[rownames(dist_phylo), ]

# Mantel tests
library(vegan)
mantel(dist_phylo, dist_geo)
mantel.partial(dist_phylo, dist_clim, dist_geo)
mantel.partial(dist_phylo[rownames(dist_soil), ], dist_soil,
  dist_geo[rownames(dist_soil), ], na.rm = T)


### end
