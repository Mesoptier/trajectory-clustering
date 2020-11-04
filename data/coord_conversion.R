#!/usr/bin/env Rscript

library(sp)
library(rgdal)

LongLatToUTM <- function(xy, zone, south=FALSE) {
    stopifnot(is.data.frame(xy))
    coordinates(xy) <- c("Longitude", "Latitude")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
    conv_str <- paste("+proj=utm +zone=", zone, if (south) " +south" else "",
        " +datum=WGS84 +units=km +no_defs", sep='')
    res <- as.data.frame(spTransform(xy, CRS(conv_str)))
    names(res) <- c("x", "y")
    return(res)
}

ConvertDir <- function(dir, pigeons) {
    dir.create(file.path(dir, "utm"), showWarnings=FALSE)
    out_ds <- file(file.path(dir, "utm", "dataset.txt"), open="wt")

    for (i in 1:length(pigeons)) {
        pigeon <- pigeons[i]
        sub_dir <- file.path(dir, pigeon)
        conn <- file(file.path(sub_dir, "dataset.txt"), open="rt")

        lines <- readLines(conn)

        dir.create(file.path(sub_dir, "utm"), showWarnings=FALSE)
        new_dataset <- file(file.path(sub_dir, "utm", "dataset.txt"))

        writeLines(lines, new_dataset)
        writeLines(lines, out_ds)

        for (j in 1:length(lines)) {
            pigeon_file <- file.path(sub_dir, lines[j])
            pigeon_table <- read.delim(pigeon_file, header=FALSE, skip=1)
            XY <- pigeon_table[,1:2]
            XY <- XY[complete.cases(XY),]
            names(XY) <- c("Longitude", "Latitude")
            utm <- LongLatToUTM(XY, "30")
            write.table(utm, file.path(sub_dir, "utm", lines[j]),
                row.names=FALSE, sep="\t", quote=FALSE)
            write.table(utm, file.path(dir, "utm", lines[j]), row.names=FALSE,
                sep="\t", quote=FALSE)
        }

        close(new_dataset)
        close(conn)
    }
    close(out_ds)
}

pigeons <- list("a55", "brc", "c17", "c35", "p29", "p39", "p94")
dir <- "Data_for_Mann_et_al_RSBL/Bladon & Church route recapping/bladon heath"
ConvertDir(dir, pigeons)

pigeons <- list("a94", "c22", "c70", "k77", "l29", "liv", "r47", "s93")
dir <- "Data_for_Mann_et_al_RSBL/Bladon & Church route recapping/church hanborough"
ConvertDir(dir, pigeons)

pigeons <- list("H22", "H27", "H30", "H35", "H38", "H41", "H42", "H71")
dir <- "Data_for_Mann_et_al_RSBL/Horspath"
ConvertDir(dir, pigeons)

pigeons <- list("H23", "H31", "H32", "H34", "H36", "H50", "H58", "H62")
dir <- "Data_for_Mann_et_al_RSBL/Weston"
ConvertDir(dir, pigeons)
