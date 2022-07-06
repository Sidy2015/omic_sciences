## code to prepare `DATASET` dataset goes here
library(utils)
load(".RData")
dfGOBP = read.delim("data/2021_Decenber30_GOID_GOBP_SGD.txt",stringsAsFactors = F,check.names = FALSE)

GOBPgmt   = readRDS("data/2021_Decenber30_GO_BP.RDS")

usethis::use_data(GOBPgmt, overwrite = TRUE)

usethis::use_data(dfGOBP, overwrite = TRUE)
