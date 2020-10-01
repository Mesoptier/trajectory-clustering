#!/usr/bin/env Rscript

library(sp)
library(rgdal)
library(mapmisc)

base <- "movebank"
dir.create(file.path(base, "projection"), showWarnings=FALSE)
conn <- file(file.path(base, "dataset.txt"), open="r")
lines <- readLines(conn)

start <- c(0, 0)
end <- c(0, 0)
count <- length(lines)
for (i in 1:count) {
    trajectory <- read.csv(file.path(base, lines[i]), header=TRUE)
    start <- start + as.numeric(head(trajectory, 1))
    end <- end + as.numeric(tail(trajectory, 1))
}
start <- start / count
end <- end / count

tpcrs <- tpeqd(rbind(start, end))
for (i in 1:count) {
    tr <- read.csv(file.path(base, lines[i]), header=TRUE)
    tr_in <- SpatialPoints(tr, proj4string=CRS("+proj=longlat +datum=WGS84"),
        bbox=NULL)
    tr_proj <- spTransform(tr_in, tpcrs)
    write.table(tr_proj, file.path(base, "projection", lines[i]),
        row.names=FALSE, col.names=FALSE)
}
