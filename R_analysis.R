### load libraries and functions
library(ape)
library(dplyr)
library(pegas)
library(phytools)
library(plyr)
library(reshape)
source("R_code/functions.R")


### load phylogenetic trees
tree_mb <- read.tree("Bayes_tree/con_50_majrule.tre")
tree_ml <- read.tree("RAxML_tree/RAxML_bipartitions.ITS_LSU_RPB2_TEF")

# drop repeated tips
tree_mb <- drop.tip(tree_mb, c("Cadophora_antarctica_FMR16056",
  "Mycochaetophora_gentianae"))
tree_ml <- drop.tip(tree_ml, c("Cadophora_antarctica_FMR16056",
  "Mycochaetophora_gentianae"))

### load strains' data
strains <- read.csv("Strains_data/strains.csv", h = T, sep = ";", row.names = 1)
refs <- read.csv("Strains_data/references.csv", h = T, sep = ";", row.names = 1)
refs <- get_country(refs)

# calculate new variables
strains$con_ratio <- strains$con_len / strains$con_wid
strains$con_vol   <- (4 / 3) * pi * strains$con_len * strains$con_wid^2 


### plot phylogenetic trees with source information
plot_fancy_tree(tree_mb, file = "Bayes_tree.pdf")
plot_fancy_tree(tree_ml, file = "RAxML_tree.pdf")

# plot map for legend
map_legend()


### compare phylogeny with OTU clustering
# subset tree to keep only isolates
tree_strains <- drop.tip(tree_mb, grep("^P[0-9]", tree_mb$tip.label,
  invert = T))
tree_strains <- ladderize(tree_strains, right = F)

# load OTU clusters at different percent identities
otus_97 <- read.csv("OTU_clustering/otu_list_97.csv", h = F, sep = "\t",
  row.names = 2)
otus_97 <- otus_97[tree_strains$tip.label, ]

otus_98 <- read.csv("OTU_clustering/otu_list_98.csv", h = F, sep = "\t",
  row.names = 2)
otus_98 <- otus_98[tree_strains$tip.label, ]

otus_99 <- read.csv("OTU_clustering/otu_list_99.csv", h = F, sep = "\t",
  row.names = 2)
otus_99 <- otus_99[tree_strains$tip.label, ]

# plot trees with OTU clusters
pdf("tree_OTUs.pdf", h = 8, w = 10)
par(mfrow = c(1, 3), xpd = T)
plot(tree_strains, show.tip.label = F)
tiplabels(otus_97, frame = "n", cex = 0.75, adj = 0, offset = 0.001,
  col = color(length(unique(otus_97)))[otus_97])
mtext("97% similarity", side = 3, adj = 0, line = -1)
plot(tree_strains, show.tip.label = F)
tiplabels(otus_98, frame = "n", cex = 0.75, adj = 0, offset = 0.001,
  col = color(length(unique(otus_98)))[otus_98])
mtext("98% similarity", side = 3, adj = 0, line = -1)
plot(tree_strains, show.tip.label = F)
tiplabels(otus_99, frame = "n", cex = 0.75, adj = 0, offset = 0.001,
  col = color(length(unique(otus_99)))[otus_99])
mtext("99% similarity", side = 3, adj = 0, line = -1)
dev.off()


### analysis of macromorphological traits
pcam <- trait_pca(file = "PCA_trait.pdf", include_qual = F)
pcam2 <- trait_pca(file = "PCA_trait_cat.pdf", include_qual = T)

trait_pca_media(file = "PCA_trait_media.pdf", include_qual = F)


pcam$ind$coord
tree <- tree_mb
tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in%
  rownames(pcam$ind$coord))])

trait_dist <- as.matrix(dist(pcam$ind$coord, diag = T, upper = T))
phylo_dist <- cophenetic(tree)
phylo_dist <- phylo_dist[rownames(trait_dist), colnames(trait_dist)]

plot(trait_dist, phylo_dist)
require(vegan)
mantel(trait_dist, phylo_dist)


### ancestral character reconstruction

# plot trees
plot_tree_cont(sqrt(strains$ncon), file = "acr_ncon.pdf")
plot_tree_cont(strains$con_len, file = "acr_con_len.pdf")
plot_tree_cont(strains$con_wid, file = "acr_con_wid.pdf")
plot_tree_cont(sqrt(strains$con_vol), file = "acr_con_vol.pdf")
plot_tree_cont(strains$con_ratio, file = "acr_con_ratio.pdf")
plot_tree_cont(strains$cma_diam, file = "acr_cma_diam.pdf")
plot_tree_cont(strains$pda_diam, file = "acr_pda_diam.pdf")
plot_tree_cont(strains$mea_diam, file = "acr_mea_diam.pdf")
plot_tree_cont(strains$CMA_pigm, file = "acr_cma_pigm.pdf")
plot_tree_cont(strains$PDA_pigm, file = "acr_pda_pigm.pdf")
plot_tree_cont(strains$MEA_pigm, file = "acr_mea_pigm.pdf")




# simulate discrete trait
disc <- strains$MEA_simp_color
names(disc) <- rownames(strains)
disc <- na.omit(disc)

tree <- drop.tip(tree_mb,
  tree_mb$tip.label[!(tree_mb$tip.label %in% names(disc))])
tree <- ladderize(tree, right = F)

# recontruct ancestral state of discrete data
cols <- setNames(palette()[1:length(unique(disc))], sort(unique(disc)))

tree$edge.length <- tree$edge.length + 0.0001
disc_fit <- ace(disc, tree, model = "ER", type = "discrete")

plotTree(tree, fsize = 0.8, ftype = "i")
nodelabels(node = 1:tree$Nnode + Ntip(tree),
  pie = disc_fit$lik.anc, piecol = cols, cex = 0.5)
tiplabels(pie = to.matrix(disc, sort(unique(disc))),
  piecol = cols, cex = 0.3)
add.simmap.legend(colors = cols, prompt = F, x = 0.9 * par()$usr[1],
  y = -max(nodeHeights(tree)), fsize = 0.8)



### ancetral state reconstruction with Ondrej's data
species_ecology <- read.csv("Strains_data/species_ecology.csv", h = T, sep = ";")
tree_mb$tip.label[tree_mb$tip.label %in% species_ecology$strain]
species_ecology$strain[!(species_ecology$strain %in% tree_mb$tip.label)]
tree_species <- drop.tip(tree_mb,
  tree_mb$tip.label[!(tree_mb$tip.label %in% species_ecology$strain)])

species <- strains$species
names(species) <- rownames(strains)
species <- species[which(names(species) != "P1615")]
species <- species[which(names(species) != "P1909")]
species <- as.data.frame(species)
species$species <- gsub(" ", "_", species$species)

refs_species <- refs$organism
names(refs_species) <- rownames(refs)
refs_species <- as.data.frame(refs_species)
colnames(refs_species) <- "species"
species <- rbind(species, refs_species)
rm(refs_species)

species <- species[tree_species$tip.label, ]
tree_species$tip.label <- species
rm(species)

tree_species <- drop.tip(tree_species,
  tree_species$tip.label[!(tree_species$tip.label %in%
  species_ecology$species)])
species_ecology <- species_ecology[species_ecology$species %in%
  tree_species$tip.label, ]
species_ecology <- droplevels(species_ecology)
rownames(species_ecology) <- species_ecology$species
species_ecology <- species_ecology[tree_species$tip.label, ]

acs_discrete <- function(disc, tree, title = NULL, cols = NULL, um = F) {
  if(is.null(cols)) {
    cols <- setNames(palette()[1:length(unique(disc))], sort(unique(disc)))
  } else {
    cols <- setNames(cols, sort(unique(disc)))
  }
  if(um == T) tree <- force.ultrametric(tree)
  tree$edge.length <- tree$edge.length + 0.0001
  tree <- ladderize(tree, right = F)
  disc_fit <- ace(disc, tree, model = "ER", type = "discrete")

  plotTree(tree, fsize = 0.8, ftype = "i", no.marging = F, offset = 0.5)
  par(fg = "transparent")
  nodelabels(node = 1:tree$Nnode + Ntip(tree),
    pie = disc_fit$lik.anc, piecol = cols, cex = 0.5)
  tiplabels(pie = to.matrix(disc, sort(unique(disc))),
    piecol = cols, cex = 0.5)
  par(fg = "black")
  legend("topleft", legend = names(cols), pch = 20, col = cols, bty = "n",
    inset = c(0, 0.05), pt.cex = 2)
  if(!is.null(title)) mtext(title, side = 3, line = -1, adj = 0, cex = 0.8)
  return(disc_fit)
}


pdf("ASR_reduced.pdf", w = 8, h = 8)
par(mfrow = c(2, 2))
acs_discrete(species_ecology$conidiogenesis, tree_species, "Conidiogenesis",
  cols = c("chocolate1", "cadetblue3", gray(0.5)), um = T)
acs_discrete(species_ecology$conidia, tree_species, "Conidia",
  cols = c("darkolivegreen3", "brown2", "darkcyan", gray(0.5)), um = T)
acs_discrete(species_ecology$conidia_type, tree_species, "Conidia type",
  cols = c("chartreuse2", "burlywood2", "cadetblue2", "gold3", gray(0.5)), um = T)
acs_discrete(species_ecology$substrate, tree_species, "Substrate",
  cols = c("chartreuse1", "brown3", "darkgoldenrod1", "black",
  "chartreuse4", gray(0.5), "darkgoldenrod4"), um = T)
dev.off()


acs_discrete(species_ecology$conidiogenesis, tree_species, "Conidiogenesis",
  cols = c("chocolate1", "cadetblue3", gray(0.5)), um = T)$lik.anc
nodelabels(1:40)
acs_discrete(species_ecology$conidia, tree_species, "Conidia",
  cols = c("darkolivegreen3", "brown2", "darkcyan", gray(0.5)), um = T)$lik.anc
nodelabels(1:40)
acs_discrete(species_ecology$conidia_type, tree_species, "Conidia type",
  cols = c("chartreuse2", "burlywood2", "cadetblue2", "gold3", gray(0.5)),
  um = T)$lik.anc
nodelabels(1:40)
acs_discrete(species_ecology$substrate, tree_species, "Substrate",
  cols = c("chartreuse1", "brown3", "darkgoldenrod1", "black",
  "chartreuse4", gray(0.5), "darkgoldenrod4"), um = T)$lik.anc
nodelabels(1:40)
acs_discrete(species_ecology$strategy, tree_species, "Strategy",
  cols = c("darkorchid1", "brown3", "chartreuse3", "darkolivegreen3",
  gray(0.5)), um = T)$lik.anc
nodelabels(1:40)
acs_discrete(species_ecology$aquatic, tree_species, "Aquatic",
  cols = c("darkgoldenrod", "aquamarine3"), um = T)$lik.anc
nodelabels(1:40)


pdf("ASR.pdf", w = 8, h = 12)
par(mfrow = c(3, 2))
acs_discrete(species_ecology$conidiogenesis, tree_species, "Conidiogenesis",
  cols = c("chocolate1", "cadetblue3", gray(0.5)), um = T)
acs_discrete(species_ecology$conidia, tree_species, "Conidia",
  cols = c("darkolivegreen3", "brown2", "darkcyan", gray(0.5)), um = T)
acs_discrete(species_ecology$conidia_type, tree_species, "Conidia type",
  cols = c("chartreuse2", "burlywood2", "cadetblue2", "gold3", gray(0.5)),
  um = T)
acs_discrete(species_ecology$substrate, tree_species, "Substrate",
  cols = c("chartreuse1", "brown3", "darkgoldenrod1", "black",
  "chartreuse4", gray(0.5), "darkgoldenrod4"), um = T)
acs_discrete(species_ecology$strategy, tree_species, "Strategy",
  cols = c("darkorchid1", "brown3", "chartreuse3", "darkolivegreen3",
  gray(0.5)), um = T)
acs_discrete(species_ecology$aquatic, tree_species, "Aquatic",
  cols = c("darkgoldenrod", "aquamarine3"), um = T)
dev.off()


### calculate phylogenetic signal for each variable
# select unique species
tree_sp <- select_species(tree_mb, only_strains = T, use_species_names = F)
tree_sp$tip.label <- as.character(strains[tree_sp$tip.label, "species"])
plot(tree_sp)
traits_quan <- strains[, c("cma_diam", "pda_diam", "mea_diam", "ncon",
  "CMA_pigm", "PDA_pigm", "MEA_pigm", "con_vol", "con_ratio")]

traits_quan_mean <- apply(traits_quan, 2,
  function(x) tapply(x, strains$species, mean, na.rm = T))
traits_quan_sd <- apply(traits_quan, 2,
  function(x) tapply(x, strains$species, sd, na.rm = T))



phylosig_traits <- matrix(NA, ncol = 2, nrow = ncol(traits_quan_mean),
  dimnames = list(colnames(tab), c("K", "P")))
rownames(phylosig_traits) <- colnames(traits_quan_mean)
for(i in 1:ncol(traits_quan_mean)){
  x <- phylosig(tree_sp, traits_quan_mean[, i], test = T, method = "K")
  phylosig_traits[i, "K"] <- x$K
  phylosig_traits[i, "P"] <- x$P
}
rm(x)
phylosig_traits
write.table(phylosig_traits, "phylosig_traits.csv", sep = ";", col.names = NA)



### distribution maps

# import GenBank BLAST
blast <- read.csv("Distributions/blast_result.txt", h = F, sep = "\t")
colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sgi")
blast$sequence_version <- unlist(lapply(strsplit(as.character(blast$sseqid),
  "|", fixed = T), '[[', 4))

# import SRA BLAST
sra_blast <- read.csv("Distributions/sra_blast_output_filtered.csv", h = F,
  sep = "\t")
colnames(sra_blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
sra_blast$sample <- unlist(lapply(strsplit(gsub(".", "|", sra_blast$sseqid,
  fixed = T), "|", fixed = T), '[[', 3))

# remove multiple detections per fungus and SRA, to reduce memory use
nrow(sra_blast)
sra_blast <- sra_blast[rownames(unique(sra_blast[, c("qseqid", "sample")])), ]
nrow(sra_blast)

# import GenBank records
data <- read.csv("Distributions/gb_records.csv", h = T, sep = "\t")
data <- unique(data)
data <- get_country(data)
data$lat_lon <- as.character(data$lat_lon)

# manually input data on sequences
modify_lat_lon("MF979577", "51.122541 N 115.382972 W") # C. interclivum
modify_lat_lon("MF979574", "51.122541 N 115.382972 W") # C. meredithiae


##### NEEDS TO BE CHECKED!!
for(i in 1:nrow(data)) {
  if(is.na(data$lat[i]) & !is.na(data$country2[i])) {
    xy <- NA
    xy <- try(get_centroid(data$country2[i]), silent = T)
    data$lon[i] <- xy[1]
    data$lat[i] <- xy[2]
  }
}
data$lat <- as.numeric(data$lat)
data$lon <- as.numeric(data$lon)
#####


# join blast data and metadata
data <- left_join(blast, data, by = "sequence_version")
data <- droplevels(data[data$pident > 98, ])


data_sra1 <- read.csv("Distributions/RunInfo.csv", h = T, sep = ",")
colnames(data_sra1) <- tolower(colnames(data_sra1))
data_sra1$sample <- as.character(data_sra1$sample)

data_sra2 <- read.csv("Distributions/SRA_metadata_tab.csv", h = T, sep = "\t")
data_sra2 <- droplevels(data_sra2[data_sra2$lat_lon != "missing", ])
data_sra2$lat <- unlist(lapply(strsplit(as.character(data_sra2$lat_lon),
  " ", fixed = T), '[[', 1))
data_sra2$lon <- unlist(lapply(strsplit(as.character(data_sra2$lat_lon),
  " ", fixed = T), '[[', 2))
data_sra2$lat <- as.numeric(data_sra2$lat)
data_sra2$lon <- as.numeric(data_sra2$lon)
data_sra <- left_join(data_sra1, data_sra2, by = "biosample")
rm(data_sra1, data_sra2)

save.image()
#gc()


head(data_sra)
data_sra <- data_sra[order(data_sra$srastudy), ]
SRA_obj <- as.data.frame(rowSums(table(data_sra$srastudy, data_sra$biosample)))
colnames(SRA_obj) <- "n_samples"
SRA_obj$reads <- tapply(data_sra$spots, data_sra$srastudy, sum)
SRA_obj$median_length <- tapply(data_sra$avglength, data_sra$srastudy, median)
SRA_obj$platform <- unique(data_sra[ , c("srastudy", "model")])$model
SRA_obj$organism <- na.omit(unique(data_sra[ , c("srastudy", "organism")]))$organism
SRA_obj$region <- c("Japan", "USA", "Global", "Europe", "USA", "Europe",
  "Europe", "Australia", "Australia")
SRA_obj$reference <- c("unpublished",
  "https://doi.org/10.1111/mec.12821",
  "https://doi.org/doi/10.1126/science.1256688",
  "https://doi.org/10.1111/nph.14873",
  "https://doi.org/10.3389/fmicb.2019.00481",
  "https://doi.org/10.1073/pnas.1710455114",
  "https://doi.org/10.1080/17550874.2018.1504332",
  "unpublished",
  "https://doi.org/10.1088%2F1748-9326%2Faadc19")

data_sra <- droplevels(data_sra[data_sra$sample %in% sra_blast$sample, ])
sra_blast <- droplevels(sra_blast[sra_blast$sample %in% data_sra$sample, ])

data_sra <- left_join(data_sra, sra_blast, by = "sample")
#data_sra <- full_join(data_sra, sra_blast, by = "sample")
data_sra <- droplevels(data_sra[data_sra$pident > 98 & data_sra$length > 200, ])
rm(sra_blast)


if (!dir.exists("Distribution_maps")) {
  dir.create("Distribution_maps", recursive = T)
}
for(i in unique(data$qseqid)) plot_distribution(i)


strains$range_100 <- get_species_ranges(pident = 100)
strains$range_99 <- get_species_ranges(pident = 99)
strains$range_98 <- get_species_ranges(pident = 98)




species <- 1:length(unique(strains$species))
#names(species) <- unique(strains$species)[c(9, 5, 10, 17, 8, 14, 1, 7, 2, 4, 16,
#  13, 11, 3, 18, 12, 15, 6)]
names(species) <- unique(strains$species)[c(14, 12, 16, 6, 15, 9, 8, 5, 10, 1,
  13, 17, 7, 4, 11, 18, 3, 2)]


pdf("ranges.pdf", w = 7, h = 3.5, pointsize = 12)
par(mar = c(8, 5, 1, 1), xpd = T, las = 1, mgp = c(3.5, 0.5, 0), tcl = -0.25)
plot(0, 0, type = "n", xlim = c(1, length(species) + 1),
  ylim = c(0, max(strains$range_98, na.rm = T) * 1.05), axes = F,
  ylab = "Range (Km)", xlab = NA)
for(i in seq(5000, 20000, 5000)) {
  lines(c(0, length(species) + 0.5), c(i, i), lty = 2, col = gray(0.5))
}
axis(1, at = 0.5:(length(species) + 0.5), labels = NA, pos = 0, lwd = 0,
  lwd.ticks = 1)
axis(1, at = c(0, (length(species) + 0.5)), labels = NA, pos = 0, lwd = 1,
  lwd.ticks = 0)
axis(2, pos = 0)
text(species, rep(-1000, length(species)), names(species), srt = 45,
  adj = 1)
for(i in rownames(strains)) {
  x <- species[as.character(strains[i, "species"])]
  points(jitter(x, 0.75), strains[i, "range_98"], pch = 16,
    col = alpha("#fff700", 0.5), cex = 1.5)
}
for(i in rownames(strains)) {
  x <- species[as.character(strains[i, "species"])]
  points(jitter(x, 0.75), strains[i, "range_99"], pch = 16,
    col = alpha("#ff9b00", 0.5), cex = 1.5)
}
for(i in rownames(strains)) {
  x <- species[as.character(strains[i, "species"])]
  points(jitter(x, 0.75), strains[i, "range_100"], pch = 16,
    col = alpha("#ff0000", 0.5), cex = 1.5)
}
dev.off()


for(i in unique(strains$species)) plot_species_distribution(i)




### effect of ecological variables
# get phylogenetic distance
dist_phylo <- cophenetic(tree_strains)
dist_phylo <- dist_phylo[rownames(strains), rownames(strains)]
var_phylo <- cmdscale(dist_phylo, k = 6)

# get morphological distance
var_morph <- pcam$ind$coord
var_morph <- var_morph[rownames(dist_phylo), ]

# get climatic distance
ecol_var <- read.csv("Strains_data/ecological_variables.csv", h = T, sep = ";",
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
#var_soil <- na.omit(var_soil)

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







#vp <- varpart(dist_phylo, var_morph_mea, var_clim, var_geo)
vp <- varpart(var_phylo[rownames(var_soil), ],
  var_morph[rownames(var_soil), ],
  var_clim[rownames(var_soil), ],
  var_geo[rownames(var_soil), ],
  var_soil)
plot(vp)

vp_test <- testVP4(vp)
write.table(vp_test, "vp_test.csv", sep = ";", col.names = NA)
pdf("euler_vars.pdf", h = 4, w = 4)
plotVpEuler4(vp, names = c("Macromorphology", "Soil", "Climate", "Space"))
dev.off()

pdf("varparts.pdf", h = 4, w = 4)
showvarparts(4)
dev.off()







#### tree for Fig. 14
sp_types <- c("P1751", "P6045", "P2794", "P1323", "P2437", "P2973",
  "P1963", "P1854", "P1176")
tree_fig <- tree_mb
tree_fig <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp_types)])

pdf("species_tree.pdf", h = 4, w = 4)
par(xpd = T)
plot(tree_fig, type = "fan", align.tip.label = T)
tiplabels(pch = 16, cex = 1.5)
dev.off()






