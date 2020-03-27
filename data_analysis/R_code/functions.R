### calculate Ancestral State Reconstruction of discrete characters ------------
acs_discrete <- function(disc, tree, title = NULL, cols = NULL, um = F) {
  require(phytools)
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


### calculate phylogenetic signal of quantitative cheracters with Blomberg's K -
calc_phylosignal <- function(tree, file = "phylosig_traits.csv") {
  require(phytools)
  tree_sp <- select_species(tree, only_strains = T, use_species_names = F)
  tree_sp$tip.label <- as.character(strains[tree_sp$tip.label, "species"])
  traits_quan <- strains[, c("cma_diam", "pda_diam", "mea_diam", "ncon",
    "CMA_pigm", "PDA_pigm", "MEA_pigm", "con_vol", "con_ratio")]
  traits_quan_mean <- apply(traits_quan, 2,
    function(x) tapply(x, strains$species, mean, na.rm = T))
  traits_quan_sd <- apply(traits_quan, 2,
    function(x) tapply(x, strains$species, sd, na.rm = T))
  phylosig_traits <- matrix(NA, ncol = 2, nrow = ncol(traits_quan_mean),
    dimnames = list(colnames(traits_quan), c("K", "P")))
  rownames(phylosig_traits) <- colnames(traits_quan_mean)
  for(i in 1:ncol(traits_quan_mean)){
    x <- phylosig(tree_sp, traits_quan_mean[, i], test = T, method = "K")
    phylosig_traits[i, "K"] <- x$K
    phylosig_traits[i, "P"] <- x$P
  }
  write.table(phylosig_traits, file, sep = ";", col.names = NA)
}


### change tree tip labels -----------------------------------------------------
change_tip_names <- function(tree) {
  # modify strain names
  strain_names <- paste(strains$species, rownames(strains))
  names(strain_names) <- rownames(strains)
  for(i in rownames(strains)) {
    if(strains[i, "type"]) strain_names[i] <- paste(strain_names[i], "*")
  }
  
  # modify references names
  refs_names <- paste(refs$organism, refs$strain, sep = "_")
  names(refs_names) <- rownames(refs)
  for(i in rownames(refs)) {
    if(refs[i, "type"]) refs_names[i] <- paste(refs_names[i], "*")
  }

  # replace tip labels
  tree$tip.label <- c(strain_names, refs_names)[tree$tip.label]
  tree$tip.label <- as.character(tree$tip.label)
  return(tree)
}


### function to generate distinct colors ---------------------------------------
# adapted from https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color <- function(n, start = 1, random = F, last_gray = T) {
  require(RColorBrewer)
  
  # 433 colors
  # x <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  # 74 colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  x <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
    rownames(qual_col_pals)))

  result <- x[start:(start + n - 1)]
  if(random) result <- sample(x, n)
  if(last_gray) result[length(result)] <- gray(0.75)
  
  return(result)  
}


### find species range in Km ---------------------------------------------------
find_species_range <- function(x, lon = "lon", lat = "lat") {
  require(sp)
  xy <- as.matrix(x[, c(lon, lat)])
  xy <- na.omit(xy)
  if(nrow(xy) > 0) {
    xy_dist <- spDists(xy, longlat = T)
    result <- max(xy_dist)
  } else {
    result <- NA
  }
  return(result)  
}


### calculate climatic distances across sampling sites -------------------------
get_clim_dist <- function(file = "data/ecological_variables.csv",
  phylo = dist_phylo) {
  ecol_var <- read.csv(file, h = T, sep = ";", row.names = 1)
  ecol_var <- ecol_var[rownames(ecol_var) %in% rownames(phylo), ]
  var_clim <- ecol_var[, grep("^bio", colnames(ecol_var))]
  dist_clim <- dist(scale(var_clim), diag = T, upper = T)
  dist_clim <- as.matrix(dist_clim)
  dist_clim <- dist_clim[rownames(phylo), rownames(phylo)]
  return(dist_clim)
}


### get standardized country information ---------------------------------------
get_country <- function(data) {
  require(rnaturalearth)

  data <- split_lat_lon(data)
  country <- sub(":.*", "", data$country)
  
  # records with coordinates but no country field
  #sel <- which(!is.na(data$lat_lon) & is.na(data$country))
  #coord <- na.omit(data[sel, c("lon", "lat")])
  #country2 <- coords2country(coord)
  #country[sel] <- country2

  # manually change some names
  country[country == "USA"] <- "United States"
  country[country == "UK"] <- "United Kingdom"
  country[country == "Viet Nam"] <- "Vietnam"
  country[country == "Cote d'Ivoire"] <- "Côte d'Ivoire"
  country[country == "Svalbard"] <- "Norway"
  country[country == "Czech Republic"] <- "Czech Rep."
  country[country == "South Korea"] <- "Korea"
  #country[country == "French Guiana"] <- "Guyana"

  # check country names against map names
  map_countries <- ne_countries(scale = "medium", returnclass = "sf")$name
  message("Countries not matching map (removed):")
  print(unique(na.omit(country[!(country %in% map_countries)])))
  
  # add data to file
  data$country2 <- country
  return(data)
}


### calculate geographic distances among sampling sites ------------------------
get_geo_dist <- function(data = strains, phylo = dist_phylo) {
  require(sp)
  dist_geo <- spDists(as.matrix(data[, c("longitude", "latitude")]),
    longlat = T)
  rownames(dist_geo) <- rownames(data)
  colnames(dist_geo) <- rownames(data)
  dist_geo <- dist_geo[rownames(dist_geo) %in% rownames(phylo), ]
  return(dist_geo)
}


### calculate morphological distance from PCAmis data --------------------------
get_morph_dist <- function(pcam, phylo = dist_phylo) {
  dist_morph <- pcam$ind$coord
  dist_morph <- dist_morph[rownames(dist_morph) %in% rownames(phylo), ]
  return(dist_morph)
}


### calculate phylogenetic dstance from tree -----------------------------------
get_phylo_dist <- function(tree, data = strains) {
  dist_phylo <- cophenetic(tree)
  dist_phylo <- dist_phylo[rownames(data), rownames(data)]
  return(dist_phylo)
}


### calculate soil distances across sampling sites -----------------------------
get_soil_dist <- function(file = "data/ecological_variables.csv",
  phylo = dist_phylo) {
  ecol_var <- read.csv(file, h = T, sep = ";", row.names = 1)
  ecol_var <- ecol_var[rownames(ecol_var) %in% rownames(phylo), ]
  var_soil <- ecol_var[, c("pH", "N", "organic_mater", "lime", "avP", "P",
    "avK", "S", "cond", "Na", "K", "Mg", "Ca")]
  dist_soil <- dist(scale(var_soil), diag = T, upper = T)
  dist_soil <- as.matrix(dist_soil)
  dist_soil <- dist_soil[rownames(dist_phylo), rownames(dist_phylo)]
  return(dist_soil)
}


### calculate species ranges for all species -----------------------------------
get_species_ranges <- function(x = strains, y = blast_nt, z = blast_sra,
  pident = 100) {
  ranges <- rep(NA, nrow(x))
  names(ranges) <- rownames(x)
  for(i in rownames(x)) {
    temp <- y[y$qseqid == i & y$pident >= pident, c("lon", "lat")]
    temp <- rbind(temp, z[z$qseqid == i & z$pident >= pident, c("lon", "lat")])
    temp <- na.omit(temp)
    temp <- unique(temp)
    ranges[i] <- find_species_range(temp)    
  }
  return(ranges)
}


### calculate Mantel test between phylogenetic and morphological distances -----
mantel_morpho <- function(tree, pca) {
  require(vegan)
  tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in%
    rownames(pca$ind$coord))])
  trait_dist <- as.matrix(dist(pca$ind$coord, diag = T, upper = T))
  phylo_dist <- cophenetic(tree)
  phylo_dist <- phylo_dist[rownames(trait_dist), colnames(trait_dist)]
  #plot(trait_dist, phylo_dist)
  return(mantel(trait_dist, phylo_dist))
}


### plot map for tree legend ---------------------------------------------------
map_legend <- function(file = "map_legend.pdf") {
  require(maps)
  europe <- c("Albania", "Andorra", "Austria", "Belgium", "Bosnia", "Bulgaria",
   "Croatia", "Czech Republic", "Denmark", "France", "Germany", "Great Britain",
   "Greece", "Hungary", "Italy", "Kosovo", "Liechtenstein", "Luxembourg",
   "Macedonia", "Moldova", "Montenegro", "Netherlands", "Poland", "Portugal",
   "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Switzerland", "UK")
  xy <- strains[, c("collection_site", "longitude", "latitude")]
  xy <- unique(xy) 
  rownames(xy) <- xy$collection_site
  xy <- xy[sort(rownames(xy)), -1]
  xy$col  <- site_colors(rownames(xy))
  #xy_proj <- mapproject(xy$lon, xy$lat, proj = "azequidistant")
  xy_proj <- xy[, c("longitude", "latitude")]
  colnames(xy_proj) <- c("x", "y")
  xy_proj$x <- jitter(xy_proj$x, amount = 0.5)
  xy_proj$y <- jitter(xy_proj$y, amount = 0.5)

  pdf(file, w = 6, h = 6)
  map("world", regions = europe, xlim = c(-20, 35), ylim = c(35, 65),
    boundary = F, interior = F, fill = T, col = gray(0.85),
    #lty = 0, wrap = T, projection = "azequidistant", lforce = "l")
    lty = 0, wrap = T, lforce = "l")
  points(xy_proj$x, xy_proj$y, pch = 16, col = xy$col,
    cex = 2)
  dev.off()
}


### compare phylogenetic tree with OTU clusters --------------------------------
phylo_vs_OTUs <- function(tree, file = "tree_OTUs.pdf") {
  require(ape)
  # load OTU clusters at different percent identities
  otus_97 <- read.csv("data/otu_list_97.csv", h = F, sep = "\t", row.names = 2)
  otus_97 <- otus_97[tree$tip.label, ]
  otus_98 <- read.csv("data/otu_list_98.csv", h = F, sep = "\t", row.names = 2)
  otus_98 <- otus_98[tree$tip.label, ]
  otus_99 <- read.csv("data/otu_list_99.csv", h = F, sep = "\t", row.names = 2)
  otus_99 <- otus_99[tree$tip.label, ]

  # plot trees with OTU clusters
  pdf(file, h = 8, w = 10)
  par(mfrow = c(1, 3), xpd = T)
  plot(tree, show.tip.label = F)
  tiplabels(otus_97, frame = "n", cex = 0.75, adj = 0, offset = 0.001,
    col = color(length(unique(otus_97)))[otus_97])
  mtext("97% similarity", side = 3, adj = 0, line = -1)
  plot(tree, show.tip.label = F)
  tiplabels(otus_98, frame = "n", cex = 0.75, adj = 0, offset = 0.001,
    col = color(length(unique(otus_98)))[otus_98])
  mtext("98% similarity", side = 3, adj = 0, line = -1)
  plot(tree, show.tip.label = F)
  tiplabels(otus_99, frame = "n", cex = 0.75, adj = 0, offset = 0.001,
    col = color(length(unique(otus_99)))[otus_99])
  mtext("99% similarity", side = 3, adj = 0, line = -1)
  dev.off()
}


### plot trees with ASR analysis of discrete characters ------------------------
plot_ASR <- function(file ="ASR.pdf", file_reduced = "ASR_reduced.pdf") {
  # all characters
  pdf(file, w = 8, h = 12)
  par(mfrow = c(3, 2))
  write.table(acs_discrete(species_ecology$conidiogenesis, tree_species,
    "Conidiogenesis", cols = c("chocolate1", "cadetblue3", gray(0.5)),
     um = T)$lik.anc, "output/ASR_conidiogenesis.csv", sep = ";",
     col.names = NA)
  write.table(acs_discrete(species_ecology$conidia, tree_species, "Conidia",
    cols = c("darkolivegreen3", "brown2", "darkcyan", gray(0.5)),
    um = T)$lik.anc, "output/ASR_conidia.csv", sep = ";", col.names = NA)
  write.table(acs_discrete(species_ecology$conidia_type, tree_species, "Conidia type",
    cols = c("chartreuse2", "burlywood2", "cadetblue2", "gold3", gray(0.5)),
    um = T)$lik.anc, "output/ASR_conidia_type.csv", sep = ";", col.names = NA)
  write.table(acs_discrete(species_ecology$substrate, tree_species, "Substrate",
    cols = c("chartreuse1", "brown3", "darkgoldenrod1", "black", "chartreuse4",
    gray(0.5), "darkgoldenrod4"), um = T)$lik.anc, "output/ASR_substrate.csv",
    sep = ";", col.names = NA)
  write.table(acs_discrete(species_ecology$strategy, tree_species, "Strategy",
    cols = c("darkorchid1", "brown3", "chartreuse3", "darkolivegreen3",
    gray(0.5)), um = T)$lik.anc, "output/ASR_life_strategy.csv", sep = ";",
    col.names = NA)
  write.table(acs_discrete(species_ecology$aquatic, tree_species, "Aquatic",
    cols = c("darkgoldenrod", "aquamarine3"), um = T)$lik.anc,
    "output/ASR_aquatic.csv", sep = ";", col.names = NA)
  dev.off()

  # selection of characters
  pdf(file_reduced, w = 8, h = 8)
  par(mfrow = c(2, 2))
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
  dev.off()
}


### plot distribution maps -----------------------------------------------------
plot_distribution <- function(f, nt = blast_nt, sra = blast_sra) {
  require(ggplot2)
  require(ggthemes)

  x <- nt[nt$qseqid == f, c("pident", "lon", "lat")]
  x <- na.omit(x)
  x$point <- rep("nt", nrow(x))
  y <- sra[sra$qseqid == f, c("pident", "lon", "lat")]
  y <- na.omit(y)
  y$point <- rep("SRA", nrow(y))
  x <- rbind(x, y)
  x <- x[order(x[, 1]),]
  theme_set(theme_map())
  world <- ne_countries(scale = "medium", returnclass = "sf")
  ggplot(data = world) + geom_sf(color = gray(0.6), fill = gray(0.9)) +
    geom_point(data = x, aes(x = lon, y = lat, shape = as.character(point),
    color = pident), size = 2.5) + scale_colour_gradient(high = "red",
    low = "yellow", breaks = c(98, 98.5, 99, 99.5, 100)) +
    labs(title = f, shape = "Source", color = "Identity (%)")
  ggsave(paste0("output/", f, "_distribution.pdf"), w = 6, h = 4)
}


### plot phylogenetic tree with isolation sites --------------------------------
plot_fancy_tree <- function(tree, file = "phylo_tree.pdf") {
  require(ape)
  pdf(file, w = 7, h = 12.5, pointsize = 12)
  tree <- ladderize(tree, right = F)
  tree_new <- change_tip_names(tree)
  tip_font <- rep(3, length(tree$tip.label))
  tip_font[grep("\\*$", change_tip_names(tree)$tip.label)] <- 4
  tree_new$tip.label <- gsub(" \\*$", "", tree_new$tip.label)
  plot(tree_new, cex = 0.5, label.offset = 0.001, no.margin = T,
    font = tip_font, edge.width = 1.2)
  strains_pos <- match(rownames(strains), tree$tip.label)
  names(strains_pos) <- rownames(strains)
  tiplabels(tip = strains_pos, pch = 16, adj = 0.6,
    col = site_colors(strains[names(strains_pos), "collection_site"]),
    cex = 0.8)
  nodelabels(round(as.numeric(tree$node.label), 2), frame = "n", cex = 0.5,
    adj = c(1.25, -0.25))
  add.scale.bar(par("usr")[2] - 0.15, 0, 0.05, lwd = 2)
  dev.off()
}


### plot customized PCAmix loadings for quantitative variables -----------------
plot_quan_loadings <- function(pcam, scale = 1){
  qual <- pcam[["levels"]]$coord
  pcam <- as.data.frame(pcam[["quanti"]]$coord)
  xl <- c(min(pcam[, 1]) * 1.5, max(pcam[, 1]) * 1.5)
  yl <- c(min(pcam[, 1]) * 1.5, max(pcam[, 2]) * 1.5)
  plot(pcam[, 1], pcam[, 2], type = "n", axes = F,
    xlab = "dim 1", ylab = "dim 2",
    xlim = xl, ylim = yl)
  lines(c(par("usr")[1], par("usr")[2]), c(0, 0), lty = 2)
  lines(c(0, 0), c(par("usr")[3], par("usr")[4]), lty = 2)
  axis(1, col = "transparent",col.ticks = gray(0.25))
  axis(2, col = "transparent",col.ticks = gray(0.25))
  box()
  arr <- cbind(rep(0, nrow(pcam)), rep(0, nrow(pcam)), pcam)
  for(i in rownames(arr)) {
    arr_col <- gray(0.25)
    arrows(arr[i, 1], arr[i, 2], arr[i, 3], arr[i, 4], length = 0.05,
      col = arr_col)
    text(arr[i, 3] * 1.05, arr[i, 4] * 1.05, i, col = arr_col, cex = 0.75)
  }
  if(!is.null(qual)) {
    par(new = T)
    point_type <- rep(16, nrow(qual))
    point_type[grep("conidiogenesis=", rownames(qual))] <- 15
    point_type[grep("conidia_attachment=", rownames(qual))] <- 18
    plot(qual[, 1] * scale, qual[, 2] * scale, type = "n", axes = F, bty = "n",
      xlab = "", ylab = "")
    points(qual[, 1], qual[, 2], pch = point_type)
    text(qual[, 1], qual[, 2], rownames(qual), cex = 0.75)
    axis(side = 3)
    axis(side = 4)
  }
}


### plot distribution maps per species -----------------------------------------
plot_species_distribution <- function(sp, nt = blast_nt, sra = blast_sra) {
  require(ggthemes)

  st <- rownames(strains[strains$species == sp, ])
  x <- nt[nt$qseqid %in% st, c("pident", "lon", "lat")]
  x <- na.omit(x)
  x$point <- rep("nt", nrow(x))
  y <- sra[sra$qseqid %in% st, c("pident", "lon", "lat")]
  y <- na.omit(y)
  y$point <- rep("SRA", nrow(y))
  x <- rbind(x, y)
  x <- x[order(x[, 1]),]
  x <- unique(x)
  theme_set(theme_bw())
  world <- ne_countries(scale = "medium", returnclass = "sf")
  ggplot(data = world) + geom_sf(color = "transparent", fill = gray(0.75)) +
    geom_point(data = x, aes(x = lon, y = lat, shape = as.character(point),
    color = pident), size = 3) + scale_colour_gradient(high = "red",
    low = "yellow") + labs(title = sp, shape = "Source", color = "Identity (%)")
  ggsave(paste0("output/", sp, "_distribution.pdf"), w = 5, h = 3)
}


### plot species distribution ranges -------------------------------------------
plot_spp_ranges <- function(file = "ranges.pdf") {
  require(Hmisc)
  species <- 1:length(unique(strains$species))
  names(species) <- unique(strains$species)[c(14, 12, 16, 6, 15, 9, 8, 5, 10, 1,
    13, 17, 7, 4, 11, 18, 3, 2)]

  pdf(file, w = 7, h = 3.5, pointsize = 12)
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
}


### plot customized PCAmix square loadings -------------------------------------
plot_sq_loadings <- function(pcam, qual){
  pcam <- as.data.frame(pcam["sqload"])
  plot(pcam[, 1], pcam[, 2], type = "n", axes = F,
    xlab = "dim 1", ylab = "dim 2",
    xlim = c(0, max(pcam[, 1]) * 1.2),
    ylim = c(0, max(pcam[, 2]) * 1.2))
  #lines(c(par("usr")[1], par("usr")[2]), c(0, 0), col = gray(0.5))
  #lines(c(0, 0), c(par("usr")[3], par("usr")[4]), col = gray(0.5))
  axis(1)
  axis(2)
  box()
  arr <- cbind(rep(0, nrow(pcam)), rep(0, nrow(pcam)), pcam)
  qual <- as.data.frame(qual)
  for(i in rownames(arr)) {
    arr_col <- gray(0.2)
    if(i %in% colnames(qual)) arr_col <- "brown3"
    arrows(arr[i, 1], arr[i, 2], arr[i, 3], arr[i, 4], length = 0.05,
      col = arr_col)
    text(arr[i, 3] * 1.05, arr[i, 4] * 1.05, i, col = arr_col, cex = 0.75)
  }
}


### plot ancestral character tree with continuous variable ---------------------
plot_tree_cont <- function(var, file = "tree_cont.pdf") {
  require(phytools)
  
  cont <- var
  names(cont) <- rownames(strains)
  cont <- na.omit(cont)

  tree <- drop.tip(tree_mb,
    tree_mb$tip.label[!(tree_mb$tip.label %in% names(cont))])
  tree <- ladderize(tree, right = F)
  tree <- force.ultrametric(tree)

  # recontruct ancestral state of continuous data
  cont_fit <- fastAnc(tree, cont, vars = T, CI = T)
  obj <- contMap(tree, cont, plot = F)

  pdf(file, w = 1.5, h = 4, pointsize = 12)
  par(xpd = T)
  plot(obj, legend = 0.7 * max(nodeHeights(tree)), outline = T,
    fsize = c(0.7, 0.9), ftype = "off", lwd = 3, sig = 0)
  tiplabels(pch = 15,
    col = species_col()[strains[tree$tip.label, "species"]], offset = 0.01)
  dev.off()
}


### select unique species from tree --------------------------------------------
select_species <- function(tree, only_strains = F, use_species_names = T) {
  species <- strains$species
  names(species) <- rownames(strains)
  species <- species[which(names(species) != "P1615")]
  species <- species[which(names(species) != "P1909")]
  species <- as.data.frame(species)
  species <- unique(species)
  species$species <- gsub(" ", "_", species$species)

  if(!only_strains) {
    # modify references names
    refs_species <- refs$organism
    names(refs_species) <- rownames(refs)
    refs_species <- as.data.frame(refs_species)
    refs_species <- unique(refs_species)
    colnames(refs_species) <- "species"
    species <- rbind(species, refs_species)
    species <- unique(species)
  }
  
  tree <- drop.tip(tree,
    tree$tip.label[!(tree$tip.label %in% rownames(species))])
  species <- species[tree$tip.label, ]
  if(use_species_names) tree$tip.label <- species
  tree <- ladderize(tree, right = F)
  return(tree)
}


### site colors ----------------------------------------------------------------
site_colors <- function(sites = NULL) {
  site_color <- c("gray10", "gray70", "gray60", "gray50", "gray40", "gray30",
    "gray20", "gray80", "darkorange", "darkorange2", "orange4", "orange3",
    "orange2", "orange", "darkorange3", "darkorange4", "springgreen3",
    "springgreen4", "seagreen4", "seagreen3", "seagreen2", "springgreen",
    "darkolivegreen", "darkolivegreen4", "royalblue4", "steelblue1",
    "slateblue", "slateblue4", "slateblue1", "steelblue3", "steelblue4",
    "royalblue1", "turquoise4", "turquoise2", "skyblue4", "skyblue3", "skyblue",
    "skyblue2", "orchid1", "orchid3", "orchid4", "darkorchid4", "tomato1",
    "tomato3", "tomato4", "firebrick1", "firebrick3", "yellow1", "yellow2",
    "yellow3", "yellow4", "yellowgreen", "#ffdd55ff")
  names(site_color) <- c("BG-007", "BG-010", "BG-011", "BG-012", "BG-013",
    "BG-014", "BG-015", "BG-023", "D-100", "D-101", "D-102", "D-103", "D-104",
    "D-105", "D-11a", "D-11b", "ES-001", "ES-002", "ES-003", "ES-004", "ES-005",
    "ES-006", "ES-010", "ES-012", "F-001", "F-002", "F-004", "F-007", "F-008",
    "F-009", "F-010", "F-011", "F-013", "F-014", "F-015", "F-021", "F-023",
    "F-024", "GR-001", "GR-002", "GR-003", "GR-004", "HR-021", "HR-022",
    "HR-023", "HR-025", "HR-028", "T-024", "T-025", "T-026", "T-027", "T-028",
    "NL-1")
  if(is.null(sites)) {
    return(site_color)
  } else {
    return(site_color[as.character(sites)])
  }
}


### subset tree with one tip per species for ASR -------------------------------
species_ASR_tree <- function(tree, data) {
  tree_species <- drop.tip(tree,
    tree$tip.label[!(tree$tip.label %in% data$strain)])

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

  species <- species[tree_species$tip.label, ]
  tree_species$tip.label <- species

  tree_species <- drop.tip(tree_species,
    tree_species$tip.label[!(tree_species$tip.label %in%
      data$species)])
  data <- data[data$species %in% tree_species$tip.label, ]
  data <- droplevels(data)
  rownames(data) <- data$species
  species_ecology <<- data[tree_species$tip.label, ]
  
  return(tree_species)
}


### define species colors ------------------------------------------------------
species_col <- function(species = strains$species) {
  # define species colors
  species_col <- alpha(color(length(unique(species))), 0.5)
  species_col <- as.character(species_col)
  names(species_col) <- unique(species)
  species_col <- species_col[sort(names(species_col))]
  return(species_col)
}


### lat_lon to separate columns ------------------------------------------------
split_lat_lon <- function(data) {
  coord <- as.character(data$lat_lon)
  lat <- as.numeric(sapply(strsplit(coord, " "), "[", 1))
  y   <- sapply(strsplit(coord, " "), "[", 2)
  lon <- as.numeric(sapply(strsplit(coord, " "), "[", 3))
  x   <- sapply(strsplit(coord, " "), "[", 4)
  y[is.na(y)] <- "-"
  for(i in 1:length(y)) if(y[i] == "S") lat[i] <- -lat[i]
  x[is.na(x)] <- "-"
  for(i in 1:length(x)) if(x[i] == "W") lon[i] <- -lon[i]
  data$lon <- lon
  data$lat <- lat
  return(data)
}


### plot trait PCA -------------------------------------------------------------
trait_pca <- function(file = "PCA_trait.pdf", include_qual = T) {
  require(Hmisc)
  require(PCAmixdata)

  # select traits
  drop_strains <- c("P1443", "P1909")
  traits_quan <- strains[-which(rownames(strains) %in% drop_strains),
    c(
      "cma_diam",
      "pda_diam",
      "mea_diam",
      "ncon",
      #"CMA_color_medium",
      #"PDA_color_medium",
      #"MEA_color_medium",
      #"CMA_aut_inhibition",
      #"PDA_aut_inhibition",
      #"MEA_aut_inhibition",
      "CMA_pigm",
      "PDA_pigm",
      "MEA_pigm",
      "con_vol",
      "con_ratio"
    )]
  traits_quan$ncon <- log(traits_quan$ncon + 1)
  traits_qual <- strains[-which(rownames(strains) %in% drop_strains),
    c(
      #"CMA_simp_color",
      #"CMA_texture_top",
      #"CMA_form",
      #"CMA_margin",
      #"CMA_elevation",
      #"PDA_simp_color",
      #"PDA_texture_top",
      #"PDA_form",
      #"PDA_margin",
      #"PDA_elevation",
      #"MEA_simp_color",
      #"MEA_texture_top",
      #"MEA_form",
      #"MEA_margin",
      #"MEA_elevation",
      "conidiogenesis",
      "conidia_attachment",
      "conidia_type"
    )]
  if(!include_qual) traits_qual <- NULL

  # define species colors
  species_col <- species_col()

  qual <- traits_qual
  quan <- traits_quan
  par(mfrow = c(2, 2))
  pcam <- PCAmix(X.quanti = quan, X.quali = qual, ndim = 2, rename.level = T,
    graph = T)

  # plot PCA
  pdf(file, h = 3, w = 6, pointsize = 12)
  par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), las = 1, mgp = c(1.5, 0.5, 0),
    tcl = -0.25, xpd = T)
  plot(pcam$ind$coord, type = "n")
  lines(c(par("usr")[1], par("usr")[2]), c(0, 0), lty = 2)
  lines(c(0, 0), c(par("usr")[3], par("usr")[4]), lty = 2)
  points(pcam$ind$coord, pch = 16, cex = 1.5,
    col = species_col[as.character(strains[rownames(traits_quan), "species"])])
  #text(pcam$ind$coord[, 1], pcam$ind$coord[, 2], rownames(pcam$ind$coord),
  #  cex = 0.5)
  #if(include_qual) {
  #  plot_sq_loadings(pcam, qual)
  #} else {
    plot_quan_loadings(pcam)
  #}
  dev.off()

  return(pcam)
} 


### plot trait PCA per medium --------------------------------------------------
trait_pca_media <- function(file = "PCA_trait.pdf", include_qual = T) {
  require(Hmisc)
  require(PCAmixdata)

  # select traits
  drop_strains <- c("P1443", "P1909")
  traits_quan <- strains[-which(rownames(strains) %in% drop_strains),
    c(
      "cma_diam",
      "pda_diam",
      "mea_diam",
      "ncon",
      #"CMA_color_medium",
      #"PDA_color_medium",
      #"MEA_color_medium",
      #"CMA_aut_inhibition",
      #"PDA_aut_inhibition",
      #"MEA_aut_inhibition",
      "CMA_pigm",
      "PDA_pigm",
      "MEA_pigm",
      "con_vol",
      "con_ratio"
    )]
  traits_quan$ncon <- log(traits_quan$ncon + 1)
  traits_qual <- strains[-which(rownames(strains) %in% drop_strains),
    c(
      "CMA_simp_color",
      "CMA_texture_top",
      "CMA_form",
      "CMA_margin",
      "CMA_elevation",
      "PDA_simp_color",
      "PDA_texture_top",
      "PDA_form",
      "PDA_margin",
      "PDA_elevation",
      "MEA_simp_color",
      "MEA_texture_top",
      "MEA_form",
      "MEA_margin",
      "MEA_elevation"
    )]
  if(!include_qual) traits_qual <- NULL

  # define species colors
  species_col <- species_col()

  # plot PCAs
  pdf(file, h = 4, w = 8, pointsize = 12)
  par(mfcol = c(2, 4), mar = c(4, 4, 1, 1), las = 1, mgp = c(1.5, 0.5, 0),
    tcl = -0.25, xpd = T)
  qual <- traits_qual[, grep("MEA", colnames(traits_qual))]
  quan <- traits_quan[, -grep("CMA", colnames(traits_quan), ignore.case = T)]
  quan <- quan[, -grep("PDA", colnames(quan), ignore.case = T)]
  pcam <- PCAmix(X.quanti = quan, X.quali = qual, ndim = 2, rename.level = T,
    graph = F)
  plot(pcam$ind$coord, pch = 16, cex = 2,
    col = species_col[as.character(strains[rownames(traits_quan), "species"])],
    main = "MEA")
  #lines(c(par("usr")[1], par("usr")[2]), c(0, 0), col = gray(0.5))
  #lines(c(0, 0), c(par("usr")[3], par("usr")[4]), col = gray(0.5))
  #text(pcam$ind$coord[, 1], pcam$ind$coord[, 2], rownames(pcam$ind$coord),
  #  cex = 0.5)
  plot_sq_loadings(pcam, qual)

  qual <- traits_qual[, grep("PDA", colnames(traits_qual))]
  quan <- traits_quan[, -grep("CMA", colnames(traits_quan), ignore.case = T)]
  quan <- quan[, -grep("MEA", colnames(quan), ignore.case = T)]
  pcam <- PCAmix(X.quanti = quan, X.quali = qual, ndim = 2, rename.level = T,
    graph = F)
  plot(pcam$ind$coord, pch = 16, cex = 2,
    col = species_col[as.character(strains[rownames(traits_quan), "species"])],
    main = "PDA")
  #lines(c(par("usr")[1], par("usr")[2]), c(0, 0), col = gray(0.5))
  #lines(c(0, 0), c(par("usr")[3], par("usr")[4]), col = gray(0.5))
  plot_sq_loadings(pcam, qual)

  qual <- traits_qual[, grep("CMA", colnames(traits_qual))]
  quan <- traits_quan[, -grep("MEA", colnames(traits_quan), ignore.case = T)]
  quan <- quan[, -grep("PDA", colnames(quan), ignore.case = T)]
  pcam <- PCAmix(X.quanti = quan, X.quali = qual, ndim = 2, rename.level = T,
    graph = F)
  plot(pcam$ind$coord, pch = 16, cex = 2,
    col = species_col[as.character(strains[rownames(traits_quan), "species"])],
    main = "CMA")
  #lines(c(par("usr")[1], par("usr")[2]), c(0, 0), col = gray(0.5))
  #lines(c(0, 0), c(par("usr")[3], par("usr")[4]), col = gray(0.5))
  plot_sq_loadings(pcam, qual)

  par(xpd = T)
  plot(pcam$ind$coord, type = "n", axes = F, ylab = NA, xlab = NA)
  legend("topleft", legend = names(species_col), pch = 16,
    col = species_col, bty = "n", pt.cex = 2)
  dev.off()
} 


### end ------------------------------------------------------------------------
