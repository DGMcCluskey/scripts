lines <- readLines("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/scripts/jhonatan.neutrophil.R")

# extract package names from library() or require()
packages <- gsub(".*\\((.*)\\).*", "\\1", grep("^(library|require)\\(", lines, value = TRUE))
packages <- gsub("[\"']", "", packages)  # remove quotes if present
packages <- unique(packages)

installed <- packages %in% rownames(installed.packages())
if(any(!installed)) {
  install.packages(packages[!installed])
}

#install.packages("later")

