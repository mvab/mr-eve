library(tidyverse)

modify_node_headers_for_neo4j <- function(x, id, idname)
{
	id_col <- which(names(x) == id)
	cl <- sapply(x, class)
	for(i in 1:length(cl))
	{
		if(cl[i] == "integer")
		{
			names(x)[i] <- paste0(names(x)[i], ":INT")
		}
		if(cl[i] == "numeric")
		{
			names(x)[i] <- paste0(names(x)[i], ":FLOAT")
		}
	}
	names(x)[id_col] <- paste0(idname, "Id:ID(", idname, ")")
	return(x)
}

modify_rel_headers_for_neo4j <- function(x, id1, id1name, id2, id2name)
{
	id1_col <- which(names(x) == id1)
	id2_col <- which(names(x) == id2)
	cl <- sapply(x, class)
	for(i in 1:length(cl))
	{
		if(cl[i] == "integer")
		{
			names(x)[i] <- paste0(names(x)[i], ":INT")
		}
		if(cl[i] == "numeric")
		{
			names(x)[i] <- paste0(names(x)[i], ":FLOAT")
		}
	}
	names(x)[id1_col] <- paste0(":START_ID(", id1name, ")")
	names(x)[id2_col] <- paste0(":END_ID(", id2name, ")")
	return(x)
}

get_gene <- function(dat, chrname, posname, radius=10000)
{
	# https://www.biostars.org/p/167818/
	# source("https://bioconductor.org/biocLite.R")
	# biocLite("Homo.sapiens")
	require(Homo.sapiens)
	require(dplyr)
	
	mycoords.gr <- data_frame(chrom=dat[[chrname]], start=dat[[posname]]-10000, end=dat[[posname]]+10000) %>% 
	makeGRangesFromDataFrame

	mycoords.gr <- mergeByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), mycoords.gr)
	mycoords.gr <- data.frame(gene_id=mycoords.gr$gene_id, gr=mycoords.gr$mycoords.gr, stringsAsFactors=FALSE) %>% as_data_frame

	mycoords.gr <- inner_join(mycoords.gr, as.data.frame(org.Hs.egSYMBOL), by="gene_id")
	mycoords.gr$pos <- mycoords.gr$gr.start + 10000

	dat <- right_join(
		subset(mycoords.gr, select=c(gr.seqnames, pos, symbol, gene_id)), 
		dat, 
		by=c("gr.seqnames"=chrname, "pos"=posname))

	return(dat)
}

ucsc_get_position <- function(snp)
{
	snp <- paste(snp, collapse="', '")
	require(RMySQL)
	message("Connecting to UCSC MySQL database")
	mydb <- dbConnect(MySQL(), user="genome", dbname="hg19", host="genome-mysql.cse.ucsc.edu")

	query <- paste0(
		"SELECT * from snp144 where name in ('", snp, "');"
	)
	message(query)
	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	# dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)
	return(d)

}


## Nodes
# - Traits
# - SNPs
# - Genes

## Relationships
# - SNP-Gene
# - SNP-Trait
# - Trait-Trait

load("../results/01/outcome_nodes.rdata")
load("../results/01/exposure_dat.rdata")

exposure_nodes <- subset(exposure_dat, !duplicated(id.exposure))
table(exposure_nodes$id.exposure %in% outcome_nodes)

load("../results/01/extract.rdata")
load("../results/01/mr.rdata")

## Traits
names(nodes)[names(nodes) == "trait"] <- "name"
nodes <- modify_node_headers_for_neo4j(nodes, "id", "trait")
write.csv(nodes, file="../results/01/upload/traits.csv", row.names=FALSE, na="")


## SNPs and genes
snp_trait <- bind_rows(l) %>% as_data_frame() %>% filter(!duplicated(SNP))

snps <- ucsc_get_position(unique(snp_trait$SNP))
snps <- subset(snps, !duplicated(name)) %>% 
	dplyr::select(name=name, chr=chrom, pos=chromStart, alleles=alleles, alleleFreqs=alleleFreqs, ucsc_func=func) %>% as_data_frame
snps <- bind_rows(snps, data.frame(name=unique(subset(snp_trait, !SNP %in% snps$name)$SNP)))
snps$id <- 1:nrow(snps)
snp_trait <- inner_join(snp_trait, dplyr::select(snps, SNP=name, snp_id=id), by="SNP")


gene_snp <- get_gene(subset(snps, !is.na(chr)), "chr", "pos")
gene <- subset(gene_snp, !duplicated(gene_id) & !is.na(gene_id)) %>%
	dplyr::select(chr=gr.seqnames, name=symbol, id=gene_id)


snps <- modify_node_headers_for_neo4j(snps, "id", "snp")
write.csv(snps, file="../results/01/upload/snps.csv", row.names=FALSE, na="")


gene <- modify_node_headers_for_neo4j(gene, "id", "gene")
write.csv(gene, file="../results/01/upload/genes.csv", row.names=FALSE, na="")


## SNP-Gene

gene_snp <- dplyr::select(gene_snp, gene_id, id) %>%
	filter(!is.na(gene_id)) %>% filter(!is.na(id))

gene_snp <- modify_rel_headers_for_neo4j(gene_snp, "gene_id", "gene", "id", "snp")
write.csv(gene_snp, file="../results/01/upload/gene_snp.csv", row.names=FALSE, na="")


## SNP-trait
snp_trait <- subset(snp_trait, select=-c(mr_keep.exposure, exposure, SNP))
names(snp_trait) <- gsub(".exposure", "", names(snp_trait))
snp_trait <- modify_rel_headers_for_neo4j(snp_trait, "snp_id", "snp", "id", "trait")
write.csv(snp_trait, file="../results/01/upload/snp_trait.csv", row.names=FALSE, na="")

## trait-trait

m1 <- bind_rows(m1) %>% 
	ungroup %>% 
	dplyr::select(-exposure, -outcome) %>%
	modify_rel_headers_for_neo4j("id.exposure", "trait", "id.outcome", "trait")

m2 <- bind_rows(m2) %>% 
	ungroup %>% 
	dplyr::select(-exposure, -outcome) %>%
	modify_rel_headers_for_neo4j("id.exposure", "trait", "id.outcome", "trait")

m1$instruments <- "All"
m2$instruments <- "Steiger"
m12 <- bind_rows(m1, m2)
write.csv(m12, "../results/01/upload/trait_trait.csv", row.names=FALSE, na="")


cmd <- paste0(
"~/neo4j-community-3.2.0/bin/neo4j-admin import", 
" --database mr-eve.db", 
" --id-type string", 
" --nodes:gene ../results/01/upload/genes.csv", 
" --nodes:snp ../results/01/upload/snps.csv", 
" --nodes:trait ../results/01/upload/traits.csv", 
" --relationships:GS ../results/01/upload/gene_snp.csv", 
" --relationships:GA ../results/01/upload/snp_trait.csv", 
" --relationships:MR ../results/01/upload/trait_trait.csv"
)
system(cmd)
