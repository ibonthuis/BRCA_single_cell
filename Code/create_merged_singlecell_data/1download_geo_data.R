library('GEOquery')

getwd()
setwd("/storage/kuijjerarea/ine/breast_met/single_cell/raw")
gse <- getGEO("GSE180286", GSEMatrix = TRUE, destdir = getwd())
show(gse)

#experimentData(gse[[2]])
filePaths = getGEOSuppFiles("GSE180286")
geo_file <- list.files("GSE180286", full.names = TRUE)

gunzip(geo_file, destname = "")
gunzip("[.]unzipped")

# List the .tar files in the directory
geo_tar_files <- list.files("GSE180286", pattern = "\\.tar$", full.names = TRUE)

# Untar each file
for (file in geo_file) {
  untar(file, exdir = "GSE180286_extracted") # Specify the output directory
}

# Check the extracted files
list.files("GSE180286_extracted", full.names = TRUE)