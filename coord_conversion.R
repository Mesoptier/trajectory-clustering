library(sp)
library(rgdal)


LongLatToUTM<-function(x,y,zone){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
}

pigeons = list("r47")
dir = "data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/church hanborough/"

for (i in 1:length(pigeons)) {
    pigeon = pigeons[i]
    sub_dir = paste(dir, pigeon, sep="")
    dataset_file = paste(sub_dir, "/dataset.txt", sep="")
    conn <- file(dataset_file, open='r')

    lines <- readLines(conn)

    new_dataset <- file(
        paste(sub_dir, "/utm/dataset.txt", sep="")
    )

    writeLines(lines, new_dataset)

    for (j in 1:length(lines)) {
        pigeon_file <- paste(c(sub_dir, "/", lines[j]), collapse="")
        pigeon_table <- read.delim(pigeon_file, header=FALSE, sep="\t")[-1,]
	print(lines[j])
        X <- pigeon_table[,1]
        X <- mapply(X, FUN=as.numeric)
        Y <- pigeon_table[,2]
        Y <- mapply(Y, FUN=as.numeric)
    
        utm <- LongLatToUTM(X, Y, "30U")
        write.table(utm, paste(c(sub_dir, "/utm/", lines[j]), collapse=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = FALSE)
	write.table(utm, paste(c(dir, "/utm/", lines[j]), collapse=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = FALSE)
    }

    close(new_dataset)
    close(conn)
}
