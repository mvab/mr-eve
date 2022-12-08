# Make trait, snp and gene nodes

library(dplyr)
library(argparse)
library(GenomicRanges)
library(myvariant)
library(ieugwasr)
source("scripts/mrever_functions.R")


# create parser object
parser <- ArgumentParser()
parser$add_argument('--idlist', required=TRUE)
parser$add_argument('--snplist', required=TRUE)
parser$add_argument('--genes', required=TRUE)
parser$add_argument('--outdir', required=TRUE)
args <- parser$parse_args()
str(args)


dir.create(file.path(args[["outdir"]], "resources", "neo4j_stage"))

## Nodes

# read instrument list

rsid <- scan(args[["snplist"]], what=character())
list1 <- ieugwasr::afl2_rsid(rsid) %>% 
  dplyr::select(chr, pos=start, ref=ref, alt=alt, rsid=rsid) %>%
  dplyr::mutate(build="hg19")

rem <- rsid[!rsid %in% list1$rsid]
list2 <- myvariant::queryVariants(rem) %>%
  as_tibble() %>%
  dplyr::select(chr=dbsnp.chrom, pos=dbsnp.hg19.start, ref=dbsnp.ref, alt=dbsnp.alt, rsid=query) %>%
  dplyr::mutate(build="hg19")

list2$complete <- apply(list2, 1, function(x) sum(x != "" & !is.na(x)))
list2 <- list2 %>% 
  dplyr::arrange(dplyr::desc(complete)) %>%
  dplyr::filter(!duplicated(rsid)) %>%
  dplyr::select(-complete)

variants <- bind_rows(list1, list2)
variants$chr[is.na(variants$chr)] <- "0"
variants$pos[is.na(variants$pos)] <- 0
dim(variants)
variantsfile <- write_simple(variants, 
                             file.path(args[["outdir"]], "resources", "neo4j_stage", "variants.csv.gz"), 
                             id1="rsid",
                             id1name="variant")

# Traits
load(args[["idlist"]])
idinfo<-  idinfo %>% dplyr::filter( id %in% c('ukb-a-248',   'ukb-b-1209',  'ukb-b-17422', 'ukb-b-19953', 'ukb-b-3768'))  #### THIS IS ADHOC SUBSET
idinfo <- idinfo %>% 
  dplyr::select(-c(mr, file)) %>%
  dplyr::mutate(
    trait=gsub('"', '', trait),
    note=gsub('"', '', trait)
  )
# idinfo <- list.files(file.path(config['outdir'], 'data')) %>% 
#     {.[!. %in% idinfo$id]} %>%
#     {tibble(id=.)} %>%
#     bind_rows(idinfo, .)

traitsfile <- write_simple(idinfo,
                           file.path(args[["outdir"]], "resources", "neo4j_stage", "traits.csv.gz"),
                           id1="id",
                           id1name="ogid"
)



# Genes
load(args[['genes']])
genesgr <- subset(genesgr, !duplicated(ensembl_gene_id)) %>% as_tibble()
names(genesgr)[1] <- "chr"
genesgr$chr <- as.character(genesgr$chr)
genesfile <- write_simple(genesgr,
                          file.path(args[["outdir"]], "resources", "neo4j_stage", "genes.csv.gz"),
                          id1="ensembl_gene_id",
                          id1name="ensembl_gene_id"
)




## Relationships


# - Variant-Gene
# Add 500kb region around each gene
load(args[['genes']])
genesgr <- subset(genesgr, !duplicated(ensembl_gene_id)) %>% as_tibble()
genesgr$origstart <- genesgr$start
genesgr$start <- pmax(0, genesgr$start - 500000)
genesgr$end <- genesgr$end + 500000
genesgr$width <- genesgr$end - genesgr$start + 1
genesgr <- GenomicRanges::GRanges(genesgr)
vars <- variants %>% 
  dplyr::filter(chr != "0" & pos != 0) %>%
  {GRanges(.$chr, IRanges(start=.$pos, end=.$pos), rsid = .$rsid)}
a <- GenomicRanges::findOverlaps(vars, genesgr) %>% as_tibble()
b <- bind_cols(
  vars[a$queryHits,] %>% as_tibble() %>% dplyr::select(rsid, varpos=start),
  genesgr[a$subjectHits,] %>% as_tibble()
)
b$tss_dist <- abs(b$varpos - b$origstart)
gv <- tibble(ensembl_gene_id = b$ensembl_gene_id, variant = b$rsid, tss_dist = b$tss_dist)

gvfile <- write_simple(gv, 
                       file.path(args[["outdir"]], "resources", "neo4j_stage", "gv.csv.gz"),
                       id1="ensembl_gene_id", id1name="ensembl_gene_id", id2="variant", id2name="variant")


