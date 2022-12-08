determine_analyses <- function(id, idlist, newidlist=NULL, what="eve")
{
  idlist <- unique(idlist)
  newidlist <- unique(newidlist)
  if(!is.null(newidlist))
  {
    if(any(idlist %in% newidlist))
    {
      stop("idlist and newidlist must be distinct")
    }
    if(id %in% idlist)
    {
      stop("id must be in newidlist")
    }
    if(!id %in% newidlist)
    {
      stop("id must be in newidlist")
    }
  } else {
    if(!id %in% idlist)
    {
      stop("id must be in idlist when newidlist not provided")
    }
  }
  
  if(what == "phewas")
  {
    if(!is.null(newidlist))
    {
      idlist <- unique(c(idlist, newidlist))
    }
    idlist <- idlist[!idlist %in% id]
    param <- dplyr::bind_rows(
      dplyr::tibble(exposure=id, outcome=idlist),
      dplyr::tibble(exposure=idlist, outcome=id)
    )
  }
  
  if(what == "exposure")
  {
    if(!is.null(newidlist))
    {
      idlist <- unique(c(idlist, newidlist))
    }
    idlist <- idlist[!idlist %in% id]
    param <- dplyr::tibble(exposure=id, outcome=idlist)
  }
  
  if(what == "outcome")
  {
    if(!is.null(newidlist))
    {
      idlist <- unique(c(idlist, newidlist))
    }
    idlist <- idlist[!idlist %in% id]
    param <- dplyr::tibble(exposure=idlist, outcome=id)
  }
  
  
  if(what == "eve")
  {
    if(!is.null(newidlist))
    {
      idlist <- idlist[!idlist %in% id]
      newidlist <- newidlist[!newidlist %in% id]
      param <- dplyr::bind_rows(
        dplyr::tibble(exposure=id, outcome=idlist),
        dplyr::tibble(exposure=id, outcome=newidlist),
        dplyr::tibble(exposure=idlist, outcome=id)
      )
    } else {
      idlist <- idlist[!idlist %in% id]
      param <- dplyr::bind_rows(
        dplyr::tibble(exposure=id, outcome=idlist)
      )
    }
  }
  param$id <- paste0(param$exposure, " -> ", param$outcome)
  return(param)
}



readml <- function(filename, ao, format="TwoSampleMR")
{
  if(format == "TwoSampleMR")
  {
    require(TwoSampleMR)
  }
  require(dplyr)
  # read ml.csv.gz of the exposure trait
  a <- read.csv(
    filename, 
    header=FALSE,
    stringsAsFactors=FALSE
  ) %>% dplyr::as_tibble(.)
  names(a) <- c("id", "SNP", "chr", "pos", "other_allele", "effect_allele", "beta", "se", "pval", "eaf", "samplesize", "ncase", "proxy_rsid", "instrument")
  
  # Fill in missing info
  id <- unique(a$id)
  r <- ao[tolower(ao$id) == tolower(id), ] # find exposure in ao
  id <- r$id
  a$id <- id
  
  stopifnot(nrow(r) == 1)
  
  a$samplesize[is.na(a$samplesize)] <- r$sample_size
  a$units <- r$unit
  a$ncase[is.na(a$ncase)] <- r$ncase
  a$ncontrol <- r$ncontrol
  a$Phenotype <- r$trait
  a$units[is.na(a$units)] <- "unknown"
  
  if(format == "TwoSampleMR")    
  {
    # Convert
    # formar exposure
    if(sum(a$instrument) > 0)
    {
      exposure_dat <- suppressWarnings(TwoSampleMR::format_data(subset(a, instrument))) %>% dplyr::as_tibble(.)
    } else {
      exposure_dat <- tibble()
    }
    # format outcome
    outcome_dat <- suppressWarnings(TwoSampleMR::format_data(subset(a, !instrument), type="outcome")) %>% dplyr::as_tibble(.)
    return(list(exposure_dat = exposure_dat, outcome_dat=outcome_dat))
  } else {
    return(a)
  }
}



write_out <- function(x, basename, header=FALSE)
{
  g <- gzfile(basename, "w")
  write.table(x, g, row.names=FALSE, col.names=FALSE, na="", sep=",")
  close(g)
  if(header) write.table(x[0,], file=paste0(basename, "_header.csv"), row.names=FALSE, col.names=TRUE, sep=",")
}



#' Modify headers for neo4j
#'
#'
#' @param x <what param does>
#' @param id <what param does>
#' @param idname <what param does>
#'
#' @export
#' @return
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

#' Modify headers for neo4j
#'
#' <full description>
#'
#' @param x <what param does>
#' @param id1 <what param does>
#' @param id1name <what param does>
#' @param id2 <what param does>
#' @param id2name <what param does>
#'
#' @export
#' @return
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

#' Write out to csv.gz, split to make more manageable files
#'
#' @param obj <what param does>
#' @param splitsize <what param does>
#' @param prefix <what param does>
#' @param id1 <what param does>
#' @param id1name <what param does>
#' @param id2=NULL <what param does>
#' @param id2name=NULL <what param does>
#'
#' @export
#' @return
write_split <- function(obj, splitsize, prefix, id1, id1name, id2=NULL, id2name=NULL)
{
  splitnum <- ceiling(length(obj) / splitsize)
  splits <- split(1:length(obj), 1:splitnum)
  nsplit <- length(splits)
  filenames <- paste0(prefix, 1:nsplit, ".csv.gz")
  lapply(1:length(splits), function(x)
  {
    message(x, " of ", length(splits))
    temp <- bind_rows(obj[splits[[x]]])
    if(is.null(id2))
    {
      temp <- modify_node_headers_for_neo4j(temp, id1, id1name)
    } else {
      temp <- modify_rel_headers_for_neo4j(temp, id1, id1name, id2, id2name)
    }
    gz1 <- gzfile(filenames[x], "w")
    if(x == 1)
    {
      write.table(temp, file=gz1, row.names=FALSE, na="", sep=",")
    } else {
      write.table(temp, file=gz1, row.names=FALSE, na="", col.names=FALSE, sep=",")
    }
    close(gz1)
  })
  return(paste(filenames, collapse=","))
}

#' Wrapper to write out files
#'
#' <full description>
#'
#' @param obj <what param does>
#' @param filename <what param does>
#' @param id1 <what param does>
#' @param id1name <what param does>
#' @param id2=NULL <what param does>
#' @param id2name=NULL <what param does>
#' @param col.names Whether to write the header line
#' @param headeronly Whether to only write the header
#'
#' @export
#' @return
write_simple <- function(obj, filename, id1, id1name, id2=NULL, id2name=NULL, col.names=TRUE, headeronly=FALSE)
{
  if(is.null(id2))
  {
    temp <- modify_node_headers_for_neo4j(obj, id1, id1name)
  } else {
    temp <- modify_rel_headers_for_neo4j(obj, id1, id1name, id2, id2name)
  }
  if(headeronly)
  {
    temp <- temp[0,]
  }
  gz1 <- gzfile(filename, "w")
  write.table(temp, file=gz1, row.names=FALSE, na="", sep=",", col.names=col.names)
  close(gz1)
  return(filename)
}
