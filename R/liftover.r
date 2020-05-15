#' Create marts for b36-38
#'
#' Slow to create these so just make once and save
#'
#' @export
#' @return Saves data to data/marts.rdata
create_marts <- function()
{
	mart38 <- useDataset(mart=useMart("ENSEMBL_MART_SNP"), dataset="hsapiens_snp")
	mart37 <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
	mart36 <- useDataset(mart=useMart("ENSEMBL_MART_SNP", host="http://may2009.archive.ensembl.org"), dataset="hsapiens_snp")
	save(mart36, mart37, mart38, file="data/marts.rdata")
}


#' Create dataset of some hm3 SNPs and their build positions
#'
#' @export
#' @return saves build_ref object
create_build_reference <- function()
{
	set.seed(314159)
	hm3 <- scan("inst/extdata/hapmap3_autosome.snplist", character()) %>% sample(., 300000)
	gwasvcf::set_bcftools()
	vcf <- gwasvcf::query_gwas("~/data/gwas/ieu-a-2/ieu-a-2.vcf.gz", rsid=hm3)
	b37 <- SummarizedExperiment::rowRanges(vcf)
	b37$rsid <- names(b37)
	GenomeInfoDb::seqlevels(b37) <- paste0("chr", GenomeInfoDb::seqlevels(b37))

	ch36 <- rtracklayer::import.chain(system.file(package="EbiDataImport", "extdata", "hg19ToHg18.over.chain"))
	b36 <- rtracklayer::liftOver(x=b37, chain=ch36) %>% unlist()

	ch38 <- rtracklayer::import.chain(system.file(package="EbiDataImport", "extdata", "hg19ToHg38.over.chain"))
	b38 <- rtracklayer::liftOver(x=b37, chain=ch38) %>% unlist()
	
	b36 <- dplyr::tibble(rsid=b36$rsid, b36=GenomicRanges::start(b36))
	b37 <- dplyr::tibble(rsid=b37$rsid, b37=GenomicRanges::start(b37))
	b38 <- dplyr::tibble(rsid=b38$rsid, b38=GenomicRanges::start(b38))

	build_ref <- merge(b36, b37, by="rsid")
	build_ref <- merge(build_ref, b38, by="rsid")
	save(build_ref, file="data/build_ref.rdata")
}


#' Lookup positions for given rsids in particular build
#'
#' @param rsid rsid
#' @param build build (36, 37 or 38)
#'
#' @export
#' @return data frame
get_positions <- function(rsid, build)
{
	data(marts)
	message("looking up build ", build)
	if(build == 36)
	{
		refsnp <- "refsnp"
	} else {
		refsnp <- "snp_filter"
	}
	biomaRt::getBM(attributes = c('refsnp_id','chr_name', 'chrom_start'), filters = c(refsnp), values = rsid, mart = get(paste0("mart", build))) %>% return()
}



#' Determines which build a set of SNps are on
#'
#' @param rsid rsid
#' @param chr chr
#' @param pos pos
#' @param build Builds to try e.g. c(37,38,36)
#'
#' @export
#' @return build if detected, or dataframe of matches if not
determine_build_biomart <- function(rsid, chr, pos, build=c(37,38,36))
{
	index <- sample(1:length(rsid), min(length(rsid), 500))
	dat <- dplyr::tibble(rsid=rsid[index], chr=chr[index], pos=pos[index])
	l <- list()
	for(i in 1:length(build))
	{
		a <- get_positions(dat$rsid, build[i])
		b <- dplyr::inner_join(a, dat, by=c("refsnp_id"="rsid"))
		n <- sum(b$chr_name == b$chr & b$chrom_start == b$pos)
		l[[i]] <- dplyr::tibble(
			build=build[i],
			prop_pos = n / nrow(b),
			prop_found = nrow(b) / nrow(dat)
		)
		if(l[[i]]$prop_pos > 0.9 & l[[i]]$prop_found > 0.3)
		{
			return(build[i])
		}
	}
	l %>% dplyr::bind_rows() %>% return()
}


#' Determine build based on reference dataset
#'
#' @param rsid rsid
#' @param chr chr
#' @param pos pos
#' @param build Builds to try e.g. c(37,38,36)
#'
#' @export
#' @return build if detected, or dataframe of matches if not
determine_build <- function(rsid, chr, pos, build=c(37,38,36))
{
	index <- sample(1:length(rsid), min(length(rsid), 1000000))
	dat <- dplyr::tibble(rsid=rsid[index], chr=chr[index], pos=pos[index])
	data(build_ref)
	dat2 <- merge(dat, build_ref, by="rsid")
	if(nrow(dat2) < 20 & nrow(dat) > 100)
	{
		message("Not enough variants in reference dataset, using biomart")

		return(determine_build_biomart(rsid, chr, pos, build))
	}
	l <- list()
	for(i in 1:length(build))
	{
		n <- sum(dat2$pos == dat2[[paste0("b", build[i])]])
		l[[i]] <- dplyr::tibble(
			build=build[i],
			prop_pos = n / nrow(dat2),
			n_found = nrow(dat2)
		)
		if(l[[i]]$prop_pos > 0.9 & l[[i]]$n_found > 20)
		{
			return(build[i])
		}
	}
	l %>% dplyr::bind_rows() %>% return()
}


#' Liftover GWAS positions
#'
#' Determine GWAS build and liftover to b37
#'
#' @param rsid rsid
#' @param chr chr
#' @param pos pos
#' @param build Builds to try e.g. c(37,38,36)
#'
#' @export
#' @return NULL if already build 37 or dataframe of b37 positions
liftover_gwas <- function(rsid, chr, pos, build=c(37, 38, 36))
{
	from <- determine_build(rsid, chr, pos, build=build)
	if(is.data.frame(from))
	{
		stop("Cannot determine build")
	} else {
		if(from == 37)
		{
			message("Already build 37")
			return(NULL)
		}
	}
	message("Build: ", from)
	if(from == 36)
	{
		path <- system.file(package="EbiDataImport", "extdata", "hg18ToHg19.over.chain")
	} else if( from == 38)
	{
		path <- system.file(package="EbiDataImport", "extdata", "hg38ToHg19.over.chain")
	} else {
		stop("from must be 36 or 38")
	}

	ch <- rtracklayer::import.chain(path)

	if(!grepl("chr", chr[1]))
	{
		chr <- paste0("chr", chr)
	}
	chr[chr == "chr23"] <- "chrX"
	chr[chr == "chr24"] <- "chrY"
	chr[chr == "chr25"] <- "chrXY"
	chr[chr == "chr26"] <- "chrM"
	chr[chr == "chrMT"] <- "chrM"


	dat <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges::IRanges(start=pos, end=pos), rsid=rsid)
	d19 <- rtracklayer::liftOver(dat, ch) %>% unlist()
	d19 <- d19 %>% dplyr::as_tibble() %>% dplyr::select(rsid=rsid, chr=seqnames, pos=start)
	d19$chr <- gsub("chr", "", d19$chr)
	return(d19)
}
