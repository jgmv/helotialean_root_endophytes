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


### plot phylogenetic tree with isolation sites --------------------------------
plot_fancy_tree <- function(tree, file = "phylo_tree.pdf") {
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


### plot map for tree legend ---------------------------------------------------
map_legend <- function(file = "map_legend.pdf") {
# plot map with countries
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


### define species colors ------------------------------------------------------
species_col <- function(species = strains$species) {
  # define species colors
  species_col <- alpha(color(length(unique(species))), 0.5)
  species_col <- as.character(species_col)
  names(species_col) <- unique(species)
  species_col <- species_col[sort(names(species_col))]
  return(species_col)
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


### calculate species ranges for all species -----------------------------------
get_species_ranges <- function(x = strains, y = data, z = data_sra,
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



### plot distribution maps -----------------------------------------------------
plot_distribution <- function(f) {
  require(ggplot2)
  require(ggthemes)

  x <- data[data$qseqid == f, c("pident", "lon", "lat")]
  x <- na.omit(x)
  x$point <- rep("nt", nrow(x))
  y <- data_sra[data_sra$qseqid == f, c("pident", "lon", "lat")]
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
  ggsave(paste0("Distribution_maps/", f, "_distribution.pdf"), w = 6, h = 4)
}


### plot distribution maps per species -----------------------------------------
plot_species_distribution <- function(sp) {
  require(ggthemes)

  st <- rownames(strains[strains$species == sp, ])
  x <- data[data$qseqid %in% st, c("pident", "lon", "lat")]
  x <- na.omit(x)
  x$point <- rep("nt", nrow(x))
  y <- data_sra[data_sra$qseqid %in% st, c("pident", "lon", "lat")]
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
  ggsave(paste0("Distribution_maps/", sp, "_distribution.pdf"), w = 5, h = 3)
}


### test significance of variation partition components (3 components) ---------
testVP3 <- function(vp, cdm = NULL) {
  # retrieve tables from vp
  y  <- eval(parse(text = vp$call[2]))
  X1 <- as.matrix(eval(parse(text = vp$call[3])))
  X2 <- as.matrix(eval(parse(text = vp$call[4])))
  X3 <- as.matrix(eval(parse(text = vp$call[5])))

  # create an output table
  tab <- rbind(vp$part[[1]][1:4], vp$part[[2]][1:4], vp$part[[3]][1:4])
  tab$percVar <- tab[, "Adj.R.square"] * 100
  tab$P <- rep(NA, nrow(tab))
  #showvarparts(3)

  if(class(y) == "dist") {
    tab[7, "P"]  <- anova.cca(dbrda(y ~ X1 + X2 + X3))$Pr[1]
    tab[1, "P"]  <- anova.cca(dbrda(y ~ X1))$Pr[1]
    tab[2, "P"]  <- anova.cca(dbrda(y ~ X2))$Pr[1]
    tab[3, "P"]  <- anova.cca(dbrda(y ~ X3))$Pr[1]
    tab[4, "P"]  <- anova.cca(dbrda(y ~ X1 + X2))$Pr[1]
    tab[5, "P"]  <- anova.cca(dbrda(y ~ X1 + X3))$Pr[1]
    tab[6, "P"]  <- anova.cca(dbrda(y ~ X2 + X3))$Pr[1]
    tab[8, "P"]  <- anova.cca(dbrda(y ~ X1 + Condition(X2) +
      Condition(X3)))$Pr[1]
    tab[9, "P"]  <- anova.cca(dbrda(y ~ X2 + Condition(X1) +
      Condition(X3)))$Pr[1]
    tab[10, "P"] <- anova.cca(dbrda(y ~ X3 + Condition(X1) +
      Condition(X2)))$Pr[1]
    tab[16, "P"] <- anova.cca(dbrda(y ~ X1 + Condition(X3)))$Pr[1]
    tab[17, "P"] <- anova.cca(dbrda(y ~ X1 + Condition(X2)))$Pr[1]
    tab[18, "P"] <- anova.cca(dbrda(y ~ X2 + Condition(X3)))$Pr[1]
    tab[19, "P"] <- anova.cca(dbrda(y ~ X2 + Condition(X1)))$Pr[1]
    tab[20, "P"] <- anova.cca(dbrda(y ~ X3 + Condition(X1)))$Pr[1]
    tab[21, "P"] <- anova.cca(dbrda(y ~ X3 + Condition(X2)))$Pr[1]
  } else {
    tab[7, "P"]  <- anova.cca(rda(y ~ X1 + X2 + X3))$Pr[1]
    tab[1, "P"]  <- anova.cca(rda(y ~ X1))$Pr[1]
    tab[2, "P"]  <- anova.cca(rda(y ~ X2))$Pr[1]
    tab[3, "P"]  <- anova.cca(rda(y ~ X3))$Pr[1]
    tab[4, "P"]  <- anova.cca(rda(y ~ X1 + X2))$Pr[1]
    tab[5, "P"]  <- anova.cca(rda(y ~ X1 + X3))$Pr[1]
    tab[6, "P"]  <- anova.cca(rda(y ~ X2 + X3))$Pr[1]
    tab[8, "P"]  <- anova.cca(rda(y ~ X1 + Condition(X2) + Condition(X3)))$Pr[1]
    tab[9, "P"]  <- anova.cca(rda(y ~ X2 + Condition(X1) + Condition(X3)))$Pr[1]
    tab[10, "P"] <- anova.cca(rda(y ~ X3 + Condition(X1) + Condition(X2)))$Pr[1]
    tab[16, "P"] <- anova.cca(rda(y ~ X1 + Condition(X3)))$Pr[1]
    tab[17, "P"] <- anova.cca(rda(y ~ X1 + Condition(X2)))$Pr[1]
    tab[18, "P"] <- anova.cca(rda(y ~ X2 + Condition(X3)))$Pr[1]
    tab[19, "P"] <- anova.cca(rda(y ~ X2 + Condition(X1)))$Pr[1]
    tab[20, "P"] <- anova.cca(rda(y ~ X3 + Condition(X1)))$Pr[1]
    tab[21, "P"] <- anova.cca(rda(y ~ X3 + Condition(X2)))$Pr[1]
  }
  return(tab)
}


### test significance of variation partition components (4 components) ---------
testVP4 <- function(vp) {
  # retrieve tables from vp
  y  <- eval(parse(text = vp$call[2]))
  X1 <- as.matrix(eval(parse(text = vp$call[3])))
  X2 <- as.matrix(eval(parse(text = vp$call[4])))
  X3 <- as.matrix(eval(parse(text = vp$call[5])))
  X4 <- as.matrix(eval(parse(text = vp$call[6])))

  # create an output table
  tab <- rbind(vp$part[[1]][1:4], vp$part[[2]][1:4], vp$part[[3]][1:4])
  tab$percVar <- tab[, "Adj.R.square"] * 100
  tab$P <- rep(NA, nrow(tab))
  #showvarparts(4)
  if(class(y) == "dist") {
    tab[15, "P"]  <- anova.cca(dbrda(y ~ X1 + X2 + X3 + X4))$Pr[1]
    tab[1, "P"]   <- anova.cca(dbrda(y ~ X1))$Pr[1]
    tab[2, "P"]   <- anova.cca(dbrda(y ~ X2))$Pr[1]
    tab[3, "P"]   <- anova.cca(dbrda(y ~ X3))$Pr[1]
    tab[4, "P"]   <- anova.cca(dbrda(y ~ X3))$Pr[1]
    tab[5, "P"]   <- anova.cca(dbrda(y ~ X1 + X2))$Pr[1]
    tab[6, "P"]   <- anova.cca(dbrda(y ~ X1 + X3))$Pr[1]
    tab[7, "P"]   <- anova.cca(dbrda(y ~ X1 + X4))$Pr[1]
    tab[8, "P"]   <- anova.cca(dbrda(y ~ X2 + X3))$Pr[1]
    tab[9, "P"]   <- anova.cca(dbrda(y ~ X2 + X4))$Pr[1]
    tab[10, "P"]  <- anova.cca(dbrda(y ~ X3 + X4))$Pr[1]
    tab[11, "P"]  <- anova.cca(dbrda(y ~ X1 + X2 + X3))$Pr[1]
    tab[12, "P"]  <- anova.cca(dbrda(y ~ X1 + X2 + X4))$Pr[1]
    tab[13, "P"]  <- anova.cca(dbrda(y ~ X1 + X3 + X4))$Pr[1]
    tab[14, "P"]  <- anova.cca(dbrda(y ~ X2 + X3 + X4))$Pr[1]
    tab[16, "P"]  <- anova.cca(dbrda(y ~ X1 + Condition(X2) +
      Condition(X3) + Condition(X4)))$Pr[1]
    tab[17, "P"]  <- anova.cca(dbrda(y ~ X2 + Condition(X1) +
      Condition(X3) + Condition(X4)))$Pr[1]
    tab[18, "P"]  <- anova.cca(dbrda(y ~ X3 + Condition(X1) +
      Condition(X2) + Condition(X4)))$Pr[1]
    tab[19, "P"]  <- anova.cca(dbrda(y ~ X4 + Condition(X1) +
      Condition(X2) + Condition(X3)))$Pr[1]
    tab[32, "P"] <- anova.cca(dbrda(y ~ X1 + Condition(X2)))$Pr[1]
    tab[33, "P"] <- anova.cca(dbrda(y ~ X1 + Condition(X3)))$Pr[1]
    tab[34, "P"] <- anova.cca(dbrda(y ~ X1 + Condition(X4)))$Pr[1]
    tab[35, "P"] <- anova.cca(dbrda(y ~ X2 + Condition(X1)))$Pr[1]
    tab[36, "P"] <- anova.cca(dbrda(y ~ X2 + Condition(X3)))$Pr[1]
    tab[37, "P"] <- anova.cca(dbrda(y ~ X2 + Condition(X4)))$Pr[1]
    tab[38, "P"] <- anova.cca(dbrda(y ~ X3 + Condition(X2)))$Pr[1]
    tab[39, "P"] <- anova.cca(dbrda(y ~ X3 + Condition(X2)))$Pr[1]
    tab[40, "P"] <- anova.cca(dbrda(y ~ X3 + Condition(X4)))$Pr[1]
    tab[41, "P"] <- anova.cca(dbrda(y ~ X4 + Condition(X1)))$Pr[1]
    tab[42, "P"] <- anova.cca(dbrda(y ~ X4 + Condition(X2)))$Pr[1]
    tab[43, "P"] <- anova.cca(dbrda(y ~ X4 + Condition(X3)))$Pr[1]
  } else {
    tab[15, "P"]  <- anova.cca(rda(y ~ X1 + X2 + X3 + X4))$Pr[1]
    tab[1, "P"]   <- anova.cca(rda(y ~ X1))$Pr[1]
    tab[2, "P"]   <- anova.cca(rda(y ~ X2))$Pr[1]
    tab[3, "P"]   <- anova.cca(rda(y ~ X3))$Pr[1]
    tab[4, "P"]   <- anova.cca(rda(y ~ X3))$Pr[1]
    tab[5, "P"]   <- anova.cca(rda(y ~ X1 + X2))$Pr[1]
    tab[6, "P"]   <- anova.cca(rda(y ~ X1 + X3))$Pr[1]
    tab[7, "P"]   <- anova.cca(rda(y ~ X1 + X4))$Pr[1]
    tab[8, "P"]   <- anova.cca(rda(y ~ X2 + X3))$Pr[1]
    tab[9, "P"]   <- anova.cca(rda(y ~ X2 + X4))$Pr[1]
    tab[10, "P"]  <- anova.cca(rda(y ~ X3 + X4))$Pr[1]
    tab[11, "P"]  <- anova.cca(rda(y ~ X1 + X2 + X3))$Pr[1]
    tab[12, "P"]  <- anova.cca(rda(y ~ X1 + X2 + X4))$Pr[1]
    tab[13, "P"]  <- anova.cca(rda(y ~ X1 + X3 + X4))$Pr[1]
    tab[14, "P"]  <- anova.cca(rda(y ~ X2 + X3 + X4))$Pr[1]
    tab[16, "P"]  <- anova.cca(rda(y ~ X1 + Condition(X2) +
      Condition(X3) + Condition(X4)))$Pr[1]
    tab[17, "P"]  <- anova.cca(rda(y ~ X2 + Condition(X1) +
      Condition(X3) + Condition(X4)))$Pr[1]
    tab[18, "P"]  <- anova.cca(rda(y ~ X3 + Condition(X1) +
      Condition(X2) + Condition(X4)))$Pr[1]
    tab[19, "P"]  <- anova.cca(rda(y ~ X4 + Condition(X1) +
      Condition(X2) + Condition(X3)))$Pr[1]
    tab[32, "P"] <- anova.cca(rda(y ~ X1 + Condition(X2)))$Pr[1]
    tab[33, "P"] <- anova.cca(rda(y ~ X1 + Condition(X3)))$Pr[1]
    tab[34, "P"] <- anova.cca(rda(y ~ X1 + Condition(X4)))$Pr[1]
    tab[35, "P"] <- anova.cca(rda(y ~ X2 + Condition(X1)))$Pr[1]
    tab[36, "P"] <- anova.cca(rda(y ~ X2 + Condition(X3)))$Pr[1]
    tab[37, "P"] <- anova.cca(rda(y ~ X2 + Condition(X4)))$Pr[1]
    tab[38, "P"] <- anova.cca(rda(y ~ X3 + Condition(X2)))$Pr[1]
    tab[39, "P"] <- anova.cca(rda(y ~ X3 + Condition(X2)))$Pr[1]
    tab[40, "P"] <- anova.cca(rda(y ~ X3 + Condition(X4)))$Pr[1]
    tab[41, "P"] <- anova.cca(rda(y ~ X4 + Condition(X1)))$Pr[1]
    tab[42, "P"] <- anova.cca(rda(y ~ X4 + Condition(X2)))$Pr[1]
    tab[43, "P"] <- anova.cca(rda(y ~ X4 + Condition(X3)))$Pr[1]
  }
  return(tab)
}


### plot Euler diagrams from variation partition object ------------------------
plotVpEuler3 <- function(vp, names = c("X1", "X2", "X3"),
  col = c("brown3", "skyblue3", "orange")) {
  require(venneuler)
  sec <- vp$part$indfract[-8, 3]
  sec <- ifelse(sec < 0, 0, sec)
  res <- vp$part$indfract[8, 3]
  names(sec)<- c(names, paste(names[1], names[2], sep = "&"),
    paste(names[2], names[3], sep = "&"), paste(names[1], names[3], sep = "&"),
    paste(names[1], names[2], names[3], sep = "&"))
  vd <- venneuler(sec, residuals = res)
  plot(vd, col = col)
  mtext(paste("Residuals = ", round(res, 2), sep = ""), 1, cex = 1)
}
plotVpEuler4 <- function(vp, names = c("X1", "X2", "X3", "X4"),
  col = c("brown3", "skyblue3", "orange", "chartreuse4")) {
  require(venneuler)
  sec <- vp$part$indfract[-16, 3]
  sec <- ifelse(sec < 0, 0, sec)
  res <- vp$part$indfract[16, 3]
  names(sec) <- c(names,
    paste(names[1], names[2], sep = "&"),
    paste(names[2], names[3], sep = "&"),
    paste(names[1], names[3], sep = "&"),
    paste(names[1], names[4], sep = "&"),
    paste(names[2], names[4], sep = "&"),
    paste(names[3], names[4], sep = "&"),
    paste(names[1], names[2], names[4], sep = "&"),# order correct?
    paste(names[1], names[2], names[3], sep = "&"),
    paste(names[2], names[3], names[4], sep = "&"),
    paste(names[1], names[3], names[4], sep = "&"),
    paste(names[1], names[2], names[3], names[4], sep = "&"))
  vd <- venneuler(sec, residuals = res)
  plot(vd, col = col)
  mtext(paste("Residuals = ", round(res, 2), sep = ""), 1, cex = 1)
}


### manually modify latitude longitude entries in gb records -------------------
modify_lat_lon <- function(accession, lat_lon) {
  sel <- grep(accession, data$accession)
  data[sel, "lat_lon"] <<- rep(lat_lon, length(sel))
}

