#!/usr/bin/env Rscript

library(sp)
library(rgdal)

UTMToLongLat <- function(xy, zone, south=FALSE) {
    stopifnot(is.data.frame(xy))
    coordinates(xy) <- c("x", "y")
    src_str <- paste("+proj=utm +zone=", zone, if (south) " +south" else "",
        " +datum=WGS84 +units=km +no_defs", sep='')
    proj4string(xy) <- CRS(src_str)
    res <- as.data.frame(spTransform(xy, CRS("+proj=longlat +datum=WGS84")))
    names(res) <- c("Longitude", "Latitude")
    return(res)
}

ConvertDir <- function(dir, pigeons) {
    for (i in 1:length(pigeons)) {
        pigeon_file <- file.path(dir, pigeons[i])
        pigeon_table <- read.csv(pigeon_file)
        longlat <- UTMToLongLat(pigeon_table, "30")
        write.table(longlat, pigeon_file, row.names=FALSE, sep=",", quote=FALSE)
    }
}

pigeons <- list("a55", "brc", "c17", "c35", "p29", "p39", "p94")
dir <- "convergence/tmp/Bladon & Church route recapping/bladon heath"
ConvertDir(dir, pigeons)

pigeons <- list("a94", "c22", "c70", "k77", "l29", "liv", "r47", "s93")
dir <- "convergence/tmp/Bladon & Church route recapping/church hanborough"
ConvertDir(dir, pigeons)

pigeons <- list("H22", "H27", "H30", "H35", "H38", "H41", "H42", "H71")
dir <- "convergence/tmp/Horspath"
ConvertDir(dir, pigeons)

pigeons <- list("H23", "H31", "H32", "H34", "H36", "H50", "H58", "H62")
dir <- "convergence/tmp/Weston"
ConvertDir(dir, pigeons)
