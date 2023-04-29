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

	ch36 <- rtracklayer::import.chain(system.file(package="GwasDataImport", "extdata", "hg19ToHg18.over.chain"))
	b36 <- rtracklayer::liftOver(x=b37, chain=ch36) %>% unlist()

	ch38 <- rtracklayer::import.chain(system.file(package="GwasDataImport", "extdata", "hg19ToHg38.over.chain"))
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
#' @param build build (36, 37 default or 38)
#' @param method "opengwas" (fastest) or "biomart"
#' @param splitsize Default 50000
#'
#' @export
#' @return data frame
#' @importFrom utils data
get_positions <- function(rsid, build=37, method=c("opengwas", "biomart")[1], splitsize=50000)
{
	if(length(rsid) > 100000)
	{
		message("This could take quite some time")
	}
	n <- length(rsid)
	chunks <- ceiling(n/splitsize)
	message("Splitting into ", chunks, " chunks of size ", splitsize, " each")
	rsidl <- split(rsid, 1:chunks)

	if(method == "opengwas")
	{
		if(build != 37)
		{
			stop("Only build 37 available on opengwas")
		}

		l <- list()
		for(i in 1:length(rsidl))
		{
			message("Chunk ", i, " of ", chunks)
			l[[i]] <- ieugwasr::afl2_rsid(rsidl[[i]])
		}
		b <- dplyr::bind_rows(l) %>%
			dplyr::select(rsid, chr, pos=start)

		missing <- rsid[! rsid %in% b$rsid]
		
		if(length(missing) > 0)
		{
			message("Missing ", length(missing), " from first pass, continuing again")
			chunks <- ceiling(length(missing)/splitsize)
			missing <- split(missing, 1:chunks)

			l <- list()
			for(i in 1:length(missing))
			{
				message("Chunk ", i, " of ", chunks)
				l[[i]] <- ieugwasr::variants_rsid(missing[[i]])
			}
			b1 <- dplyr::bind_rows(l) %>%
				dplyr::select(rsid=query, chr, pos)
			b <- dplyr::bind_rows(b, b1)
		}
	} else if(method == "biomart"){
		data(marts)
		message("looking up build ", build)
		refsnp <- ifelse(build == 36, "refsnp", "snp_filter")
		l <- list()
		for(i in 1:length(rsidl))
		{
			message("Chunk ", i, " of ", chunks)
			l[[i]] <- biomaRt::getBM(attributes = c('refsnp_id','chr_name', 'chrom_start'), filters = c(refsnp), values = rsidl[[i]], mart = get(paste0("mart", build))) %>%
				dplyr::select(rsid=refsnp_id, chr=chr_name, pos=chrom_start)
		}
	} else {
		stop("Method must be 'opengwas' or 'biomart'")
	}
	m <- match(rsid, b$rsid)
	b <- b[m, ]
	message("Found ", nrow(b), " of ", length(rsid), " rsids")
	return(b)
}



#' Determines which build a set of SNPs are on
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
		b <- dplyr::inner_join(a, dat, by=rsid)
		n <- sum(b$chr.x == b$chr.y & b$pos.x == b$pos.y)
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
#' @param fallback Whether to try "position" (fast) or "biomart" (more accurate if you have rsids) based approaches instead
#'
#' @export
#' @return build if detected, or dataframe of matches if not
determine_build <- function(rsid, chr, pos, build=c(37,38,36), fallback="position")
{
	index <- sample(1:length(rsid), min(length(rsid), 1000000))
	dat <- dplyr::tibble(rsid=rsid[index], chr=chr[index], pos=pos[index])
	data(build_ref)
	dat2 <- merge(dat, build_ref, by="rsid")
	if(nrow(dat) > 100000)
	{
		if(nrow(dat2) < 20)
		{
			message("Not enough variants in reference dataset")

			if(fallback == "position")
			{
				message("Matching by position")
				return(determine_build_position(pos, build))
			}
			if(fallback == "biomart")
			{
				message("Trying biomaRt")
				return(determine_build_biomart(rsid, chr, pos, build))
			} else {
				return(data.frame())
			}
		}
	} else {
		if(nrow(dat2) < 20)
		{
			message("Failed to match in reference dataset")
			if(biomart_fallback)
			{
				message("Trying biomaRt")
				return(determine_build_biomart(rsid, chr, pos, build))
			} else {
				return(data.frame())
			}
		}
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


#' Determine build just from positions
#'
#' A bit sketchy but computationally fast - just assumes that there will be at least 50x more position matches in the true build than either of the others.
#'
#' @param pos Vector of positions
#' @param build c(37,38,36)
#' @param threshold how many times more in the true build than the others. default = 50 
#'
#' @export
#' @return build or if not determined then dataframe
determine_build_position <- function(pos, build=c(37,38,36))
{
	stopifnot(length(build) == 3)
	data(build_ref)
	l <- list()
	for(i in 1:length(build))
	{
		n <- sum(pos %in% build_ref[[paste0("b", build[i])]])
		l[[i]] <- dplyr::tibble(
			build=build[i],
			prop_pos = n / length(pos),
			n_found = n
		)
	}
	l <- dplyr::bind_rows(l)
	print(l)
	m <- which.max(l$n_found)
	return(l$build[m])
}


#' Liftover GWAS positions
#'
#' Determine GWAS build and liftover to required build
#'
#' @param dat Data frame with chr, pos, snp name, effect allele, non-effect allele columns
#' @param build The possible builds to check data against Default = c(37,38,26)
#' @param to Which build to lift over to. Default=37
#' @param chr_col Name of chromosome column name. Required
#' @param pos_col Name of position column name. Required
#' @param snp_col Name of SNP column name. Optional. Uses less certain method of matching if not available
#' @param ea_col Name of effect allele column name. Optional. Might lead to duplicated rows if not presented
#' @param oa_col Name of other allele column name. Optional. Might lead to duplicated rows if not presented
#'
#' @export
#' @return Data frame
liftover_gwas <- function(dat, build=c(37,38,36), to=37, chr_col="chr", pos_col="pos", snp_col="snp", ea_col="ea", oa_col="oa", build_fallback="position")
{
	if(is.null(snp_col))
	{
		message("Only using position")
		from <- determine_build_position(dat[[pos_col]], build=build)
	} else {
		message("Using rsid")
		from <- determine_build(dat[[snp_col]], dat[[chr_col]], dat[[pos_col]], build=build, fallback="position")
	}
	if(is.data.frame(from))
	{
		print(from)
		stop("Cannot determine build")
	} else {
		if(from == to)
		{
			message("Already build ", to)
			return(dat)
		}
	}
	message("Lifting build: ", from, " to ", to)

	tab <- dplyr::tibble(build=c(36,37,38), name=c("Hg18", "Hg19", "Hg38"))

	path <- system.file(package="GwasDataImport", "extdata", paste0(
		tolower(tab$name[tab$build==from]),
		"To",
		tab$name[tab$build==to],
		".over.chain"
	))
	stopifnot(file.exists(path))

	message("Loading chainfile")
	ch <- rtracklayer::import.chain(path)

	message("Converting chromosome codings")
	if(!grepl("chr", dat[[chr_col]][1]))
	{
		dat[[chr_col]] <- paste0("chr", dat[[chr_col]])
	}
	dat[[chr_col]][dat[[chr_col]] == "chr23"] <- "chrX"
	dat[[chr_col]][dat[[chr_col]] == "chr24"] <- "chrY"
	dat[[chr_col]][dat[[chr_col]] == "chr25"] <- "chrXY"
	dat[[chr_col]][dat[[chr_col]] == "chr26"] <- "chrM"
	dat[[chr_col]][dat[[chr_col]] == "chrMT"] <- "chrM"


	message("Organising")
	datg <- GenomicRanges::GRanges(
		seqnames=dat[[chr_col]], 
		ranges=IRanges::IRanges(start=dat[[pos_col]], end=dat[[pos_col]]), 
		ind=1:nrow(dat)
	)

	message("Lifting")
	d19 <- rtracklayer::liftOver(datg, ch) %>% unlist()
	message("Organising again")
	dat <- dat[d19$ind,]
	dat[[chr_col]] <- d19@seqnames
	dat[[pos_col]] <- d19@ranges@start
	dat[[chr_col]] <- gsub("chr", "", dat[[chr_col]])
	message("Reordering")
	dat <- dat[order(dat[[chr_col]], dat[[pos_col]]), ]

	if(!is.null(ea_col) & !is.na(oa_col))
	{
		message("Removing duplicates")
		nom <- names(dat)
		if(is.numeric(chr_col)) chr_col <- nom[chr_col]
		if(is.numeric(pos_col)) pos_col <- nom[pos_col]
		if(is.numeric(ea_col)) ea_col <- nom[ea_col]
		if(is.numeric(oa_col)) oa_col <- nom[oa_col]

		dat <- distinct(dat, .data[[chr_col]], .data[[pos_col]], .data[[ea_col]], .data[[oa_col]], .keep_all=TRUE)
	}

	message("Done")
	return(dat)
}


liftover_gwas_old <- function(dat, build=c(37,38,36), to=37, chr_col="chr", pos_col="pos", snp_col="snp", ea_col="ea", oa_col="oa", build_fallback="position")
{
	if(is.null(snp_col))
	{
		message("Only using position")
		from <- determine_build_position(dat[[pos_col]], build=build)
	} else {
		message("Using rsid")
		from <- determine_build(dat[[snp_col]], dat[[chr_col]], dat[[pos_col]], build=build, fallback="position")
	}
	if(is.data.frame(from))
	{
		print(from)
		stop("Cannot determine build")
	} else {
		if(from == to)
		{
			message("Already build ", to)
			return(dat)
		}
	}
	message("Lifting build: ", from, " to ", to)

	tab <- dplyr::tibble(build=c(36,37,38), name=c("Hg18", "Hg19", "Hg38"))

	path <- system.file(package="GwasDataImport", "extdata", paste0(
		tolower(tab$name[tab$build==from]),
		"To",
		tab$name[tab$build==to],
		".over.chain"
	))
	stopifnot(file.exists(path))

	message("Loading chainfile")
	ch <- rtracklayer::import.chain(path)

	message("Converting chromosome codings")
	if(!grepl("chr", dat[[chr_col]][1]))
	{
		dat[[chr_col]] <- paste0("chr", dat[[chr_col]])
	}
	dat[[chr_col]][dat[[chr_col]] == "chr23"] <- "chrX"
	dat[[chr_col]][dat[[chr_col]] == "chr24"] <- "chrY"
	dat[[chr_col]][dat[[chr_col]] == "chr25"] <- "chrXY"
	dat[[chr_col]][dat[[chr_col]] == "chr26"] <- "chrM"
	dat[[chr_col]][dat[[chr_col]] == "chrMT"] <- "chrM"


	message("Organising")
	datg <- GenomicRanges::GRanges(seqnames=dat[[chr_col]], ranges=IRanges::IRanges(start=dat[[pos_col]], end=dat[[pos_col]]), LIFTOVERCHRPOS=paste0(dat[[chr_col]], ":", dat[[pos_col]]))

	message("Lifting")
	d19 <- rtracklayer::liftOver(datg, ch) %>% unlist()
	message("Organising again")
	d19 <- d19 %>% dplyr::as_tibble() %>% dplyr::select(LIFTOVERCHRPOS=LIFTOVERCHRPOS, LIFTOVERCHR=seqnames, LIFTOVERPOS=start)
	dat$LIFTOVERCHRPOS <- paste0(dat[[chr_col]], ":", dat[[pos_col]])
	dat <- merge(dat, d19, by="LIFTOVERCHRPOS")
	dat[[chr_col]] <- dat$LIFTOVERCHR
	dat[[pos_col]] <- dat$LIFTOVERPOS
	dat <- subset(dat, select=-c(LIFTOVERCHRPOS, LIFTOVERCHR, LIFTOVERPOS))
	dat[[chr_col]] <- gsub("chr", "", dat[[chr_col]])
	dat <- dat[order(dat[[chr_col]], dat[[pos_col]]), ] %>% dplyr::as_tibble()

	if(!is.null(ea_col) & !is.na(oa_col))
	{
		dat$code <- paste(dat[[chr_col]], dat[[pos_col]], dat[[ea_col]], dat[[oa_col]])
		dat <- subset(dat, !duplicated(code))
		dat <- subset(dat, select=-c(code))		
	}

	message("Done")
	return(dat)
}
