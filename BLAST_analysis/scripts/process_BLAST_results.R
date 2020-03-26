#!/usr/bin/Rscript

# load required libraries and scripts 
library(dplyr)
source("scripts/R_functions.R")

# import GenBank BLAST results
blast <- read.csv("output/blast_result.txt", h = F, sep = "\t")
colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sgi")
blast$sequence_version <- unlist(lapply(strsplit(as.character(blast$sseqid),
  "|", fixed = T), '[[', 4))

# import SRA BLAST results
sra_blast <- read.csv("output/sra_blast_output_filtered.csv", h = F,
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
data <- read.csv("output/gb_records.csv", h = T, sep = "|")
data <- unique(data)
data <- get_country(data)
data$lat_lon <- as.character(data$lat_lon)

# manually input data on some sequences
modify_lat_lon("MF979577", "51.122541 N 115.382972 W") # C. interclivum
modify_lat_lon("MF979574", "51.122541 N 115.382972 W") # C. meredithiae

# add coordinate data (centroids) for items with country names only
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

# join blast data and metadata
data <- left_join(blast, data, by = "sequence_version")
data <- droplevels(data[data$pident > 98, ])

# prepare SRA data
data_sra1 <- read.csv("data/RunInfo.csv", h = T, sep = ",")
colnames(data_sra1) <- tolower(colnames(data_sra1))
data_sra1$sample <- as.character(data_sra1$sample)

data_sra2 <- read.csv("data/SRA_metadata_tab.csv", h = T, sep = "\t")
data_sra2 <- droplevels(data_sra2[data_sra2$lat_lon != "missing", ])
data_sra2$lat <- unlist(lapply(strsplit(as.character(data_sra2$lat_lon),
  " ", fixed = T), '[[', 1))
data_sra2$lon <- unlist(lapply(strsplit(as.character(data_sra2$lat_lon),
  " ", fixed = T), '[[', 2))
data_sra2$lat <- as.numeric(data_sra2$lat)
data_sra2$lon <- as.numeric(data_sra2$lon)

data_sra <- left_join(data_sra1, data_sra2, by = "biosample")
rm(data_sra1, data_sra2)
data_sra <- data_sra[order(data_sra$srastudy), ]
data_sra <- droplevels(data_sra[data_sra$sample %in% sra_blast$sample, ])
sra_blast <- droplevels(sra_blast[sra_blast$sample %in% data_sra$sample, ])
data_sra <- left_join(data_sra, sra_blast, by = "sample")
#data_sra <- full_join(data_sra, sra_blast, by = "sample")
data_sra <- droplevels(data_sra[data_sra$pident > 98 & data_sra$length > 200, ])
rm(sra_blast)

# output tables
write.table(data, "output/blast_nt.csv", sep = ";", col.names = NA)
write.table(data_sra, "output/blast_sra.csv", sep = ";", col.names = NA)

### end
