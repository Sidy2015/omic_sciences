# installing the R packages required for class after
# installing the latest versions of R, Rstudio 
# and Bioconductor

bioc = readLines("_R_package_list.txt")

BiocManager::install(bioc)

git = readLines("_R_package_github.txt")

git = paste0("genomicsclass/",git)

devtools::install_github(git)