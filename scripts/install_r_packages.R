list.of.packages<-c("remotes", "BiocManager",
		    "dplyr", "argparse", "parallel", "randomForest", 
                    'BiocParallel', "GenomicRanges", "biomaRt", "myvariant",
                    "ieugwasr", "genetics.binaRies", "gwasvcf", "TwoSampleMR")
missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(paste0("===== Packages to install: ", missing.packages, " ====="))


# packages for installing packages
cran_mirror = "https://www.stats.bris.ac.uk/R"
if ("remotes" %in% missing.packages){install.packages("remotes", repos = cran_mirror)}
if ("BiocManager" %in% missing.packages){install.packages("BiocManager", repos = cran_mirror)}

# CRAN packages
if ("dplyr" %in% missing.packages){install.packages("dplyr", repos = cran_mirror)}
if ("argparse" %in% missing.packages){install.packages("argparse", repos = cran_mirror)}
if ("parallel" %in% missing.packages){install.packages("parallel", repos = cran_mirror)}
if ("randomForest" %in% missing.packages){install.packages("randomForest", repos = cran_mirror)}

# bioconductor packagaes
if ("BiocParallel" %in% missing.packages){BiocManager::install("BiocParallel") }# make sure this gets installed - does not work with R 4.2
if ("GenomicRanges" %in% missing.packages){BiocManager::install("GenomicRanges") }
if ("biomaRt" %in% missing.packages){BiocManager::install("biomaRt")}
if ("myvariant" %in% missing.packages){BiocManager::install("myvariant")}

# github packages
if ("ieugwasr" %in% missing.packages){remotes::install_github("mrcieu/ieugwasr")}
if ("genetics.binaRies" %in% missing.packages){remotes::install_github("MRCIEU/genetics.binaRies")}
if ("gwasvcf" %in% missing.packages){remotes::install_github("mrcieu/gwasvcf", upgrade=F) }# don't try to update dependencies
if ("TwoSampleMR" %in% missing.packages){remotes::install_github("MRCIEU/TwoSampleMR")}

missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(paste0("===== Not installed: ", missing.packages, " ====="))

if (length(missing.packages)==0){
  df <- data.frame(list.of.packages)
  write.csv(df, file = "/user/home/ny19205/scratch/mr-eve/r_packages_installed.csv", quote = F, row.names = F)
}

