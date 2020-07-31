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

pigeons <- list("a55", "brc", "c17", "c35", "p29", "p39", "p94")
dir <- "data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath"
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
        write.table(utm, file.path(sub_dir, "utm", lines[j]), row.names=FALSE,
            sep="\t", quote=FALSE)
        write.table(utm, file.path(dir, "utm", lines[j]), row.names=FALSE,
            sep="\t", quote=FALSE)
    }

    close(new_dataset)
    close(conn)
}
close(out_ds)
