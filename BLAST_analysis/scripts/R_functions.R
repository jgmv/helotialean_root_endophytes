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


### manually modify latitude longitude entries in gb records -------------------
modify_lat_lon <- function(accession, lat_lon) {
  sel <- grep(accession, data$accession)
  data[sel, "lat_lon"] <<- rep(lat_lon, length(sel))
}


### end ------------------------------------------------------------------------
