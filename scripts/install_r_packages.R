#install.packages("remotes", repos = "https://www.stats.bris.ac.uk/R")
#install.packages("BiocManager", repos = "https://www.stats.bris.ac.uk/R")

#install.packages("dplyr", repos = "https://www.stats.bris.ac.uk/R")
#install.packages("argparse", repos = "https://www.stats.bris.ac.uk/R")

#BiocManager::install("BiocParallel") # make sure this gets installed - does not work with R 4.2
#BiocManager::install("GenomicRanges")
#BiocManager::install("biomaRt")


#remotes::install_github("mrcieu/ieugwasr")
#remotes::install_github("MRCIEU/genetics.binaRies")


#remotes::install_github("mrcieu/gwasvcf", upgrade=F) # don't try to update dependencies
#remotes::install_github("MRCIEU/TwoSampleMR")

#install.packages("parallel", repos = "https://www.stats.bris.ac.uk/R")
#install.packages("randomForest", repos = "https://www.stats.bris.ac.uk/R")

BiocManager::install("myvariant")

list.of.packages<-c("dplyr", "argparse", "parallel", "randomForest", "GenomicRanges", "biomaRt", "myvariant",  "ieugwasr", "genetics.binaRies", "gwasvcf", "TwoSampleMR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(paste0("===== Not installed: ", new.packages, " ====="))
