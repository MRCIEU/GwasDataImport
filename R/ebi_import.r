#' Obtain EBI Dataset class
#'
#' @export
EbiDataset <- R6::R6Class("EbiDataset", list(
	inherit = Dataset,

	ebi_id = NULL,
	traitname = NULL,
	ftp_path = NULL,
	or_flag = NULL,

	#' @description
	#' Initialise
	#' @param ebi_id e.g. GCST004426
	#' @param wd=tempdir() <what param does>
	#' @param ftp_path=NULL <what param does>
	#' @param igd_id Defaults to "ebi-a-<ebi_id>"
	#'
	#' @return new ObtainEbiDataset object
	initialize = function(ebi_id, wd=tempdir(), ftp_path=NULL, igd_id=paste0("ebi-a-", ebi_id), traitname=NULL)
	{
		self$ebi_id <- ebi_id
		self$igd_id <- igd_id
		self$check_id(igd_id)
		self$set_wd(wd)
		self$ftp_path <- ftp_path
		self$traitname <- traitname
	},


	#' @description
	#' Download
	#' @param ftp_path=self$ftp_path <what param does>
	#' @param ftp_url=options()$ebi_ftp_url <what param does>
	#' @param outdir=self$wd <what param does>
	download_dataset = function(ftp_path=self$ftp_path, ftp_url=options()$ebi_ftp_url, outdir=self$wd)
	{
		dir.create(self$wd, recursive=TRUE, showWarnings=FALSE)
		b <- basename(ftp_path)
		filename <- file.path(outdir, b)
		ftp <- file.path(ftp_url, ftp_path)
		cmd <- glue::glue("wget -O {filename} {ftp}")
		system(cmd)
		self$filename <- filename
	},

	#' @description
	#' format dataset
	#'
	#' @param filename=self$filename <what param does>
	format_dataset = function(filename=self$filename)
	{
		keep_cols <- c("hm_rsid", "hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value")

		a <- data.table::fread(filename, header=TRUE)
		stopifnot(all(keep_cols %in% names(a)))

		or_flag <- FALSE
		if(any(!is.na(a[["hm_odds_ratio"]])))
		{
			or_flag <- TRUE
			if(! any(!is.na(a[["hm_beta"]]) & !is.na(a[["standard_error"]])))
			{
				a[["hm_beta"]] <- log(a[["hm_odds_ratio"]])
				a[["standard_error"]] <- (log(a[["hm_ci_upper"]]) - log(a[["hm_odds_ratio"]])) / 1.96
			}
		}

		out <- a %>% 
			dplyr::select(tidyselect::all_of(keep_cols)) %>%
			subset(., !is.na(hm_chrom) & !is.na(hm_pos) & !is.na(hm_effect_allele) & !is.na(hm_beta) & !is.na(standard_error))

		stopifnot(nrow(out) > 0)

		message("Determining build")
		names(out)[names(out) == "hm_rsid"] <- "rsid"
		names(out)[names(out) == "hm_chrom"] <- "chr"
		names(out)[names(out) == "hm_pos"] <- "pos"
		out <- liftover_gwas(out)
		if(!is.null(d19))
		{
			out <- merge(out, d19, by.x="hm_rsid", by.y="rsid")
			out$hm_chrom <- out$chr
			out$hm_pos <- out$pos
			out <- dplyr::select(out, tidyselect::all_of(keep_cols))
		}
		stopifnot(nrow(out) > 0)

		gwas_out <- paste0(filename, ".format.gz")
		zz <- gzfile(gwas_out, "w")
		write.table(out, file=zz, row=F, col=TRUE, qu=FALSE, na="")
		close(zz)

		self$nsnp_read <- nrow(a)
		self$nsnp <- nrow(out)
		self$gwas_out <- gwas_out
		self$or_flag <- or_flag
	},

	#' @description
	#' Download and parse metadata
	#' @param ebi_id=self$ebi_id <what param does>
	#' @param or_flag=self$or_flag <what param does>
	#' @param igd_id=NULL <what param does>
	#' @param units=NULL <what param does>
	#' @param sex="NA" <what param does>
	#' @param category="NA" <what param does>
	#' @param subcategory="NA" <what param does>
	#' @param build="HG19/GRCh37" <what param does>
	#' @param group_name="public" <what param does>
	organise_metadata = function(ebi_id=self$ebi_id, or_flag=self$or_flag, igd_id=self$igd_id, units=NULL, sex="NA", category="NA", subcategory="NA", build="HG19/GRCh37", group_name="public", traitname=self$traitname)
	{
		l <- list()
		j <- jsonlite::read_json(paste0(options()$ebi_api, ebi_id))

		l[["id"]] <- igd_id
		if(is.null(traitname))
		{
			l[["trait"]] <- j[["diseaseTrait"]][["trait"]]
		} else {
			l[["trait"]] <- traitname
		}
		l[["note"]] <- ""
		if(or_flag) l[["note"]] <- paste0(l[["note"]], "beta+se converted from OR+CI; ")
		if(!is.null(j[["studyDesignComment"]])) paste0(l[["note"]], j[["studyDesignComment"]], "; ")
		l[["pmid"]] <- j[["publicationInfo"]][["pubmedId"]]
		l[["year"]] <- j[["publicationInfo"]][["publicationDate"]]
		if(!is.null(l[["year"]])) l[["year"]] <- strsplit(l[["year"]], split="-")[[1]][1]
		l[["author"]] <- j[["publicationInfo"]][["author"]][["fullname"]]

		anc <- j[["ancestries"]]
		g <- j[["ancestries"]][[which(sapply(anc, function(x) x$type == "initial"))[1]]][["ancestralGroups"]]
		if(length(g) == 1)
		{
			l[["population"]] <- g[[1]]$ancestralGroup
		} else {
			l[["population"]] <- "Mixed"
		}

		l[["sample_size"]] <- j[["ancestries"]][[which(sapply(anc, function(x) x$type == "initial"))[1]]][["numberOfIndividuals"]]

		if(grepl("cases", j[["initialSampleSize"]]))
		{
			n <- j[["initialSampleSize"]] %>%
			gsub(",", "", .) %>%
			gsub("Up to ", "", .) %>%
			strsplit(., " ") %>%
			unlist() %>%
			as.numeric() %>%
			{.[!is.na(.)]}
			l[["ncase"]] <- n[1]
			l[["ncontrol"]] <- n[2]
			l[["unit"]] <- "logOR"
		} else {
			l[["ncase"]] <- NA
			l[["ncontrol"]] <- NA
			l[["unit"]] <- NA
		}

		if(!is.null(units))
		{
			l[["unit"]] <- units
		}

		l[["or_flag"]] <- or_flag
		l[["build"]] <- build
		l[["group_name"]] <- group_name
		l[["nsnp_stated"]] <- j[["snpCount"]]
		l[["nsnp_read"]] <- self$nsnp_read
		l[["nsnp"]] <- self$nsnp
		l[["efo"]] <- httr::GET(j[["_links"]][["efoTraits"]][["href"]]) %>% httr::content(., encoding="text") %>% 
			{.[["_embedded"]]} %>% 
			{.[["efoTraits"]]} %>%
			sapply(., function(x) x[["shortForm"]]) %>% paste(., collapse=",")
		l[["sex"]] <- sex
		l[["mr"]] <- 1
		l[["category"]] <- category
		l[["subcategory"]] <- subcategory
		self$metadata <- l

		m <- list()
		m[["id"]] <- l[["id"]]
		m[["gwas_file"]] <- basename(self$gwas_out)
		m[["snp_col"]] <- 1
		m[["chr_col"]] <- 2
		m[["pos_col"]] <- 3
		m[["ea_col"]] <- 4
		m[["oa_col"]] <- 5
		m[["eaf_col"]] <- 6
		m[["beta_col"]] <- 7
		m[["se_col"]] <- 8
		m[["pval_col"]] <- 9
		m[["header"]] <- "True"
		m[["gzipped"]] <- "True"
		m[["delimiter"]] <- "space"
		m[["cohort_cases"]] <- l[["ncase"]]
		m[["cohort_controls"]] <- ifelse(is.na(l[["ncase"]]), l[["sample_size"]], l[["ncontrol"]])
		self$datainfo <- m
	},

	pipeline = function()
	{
		message("Downloading")		
		o <- try(self$download_dataset())
		if('try-error' %in% class(o))
		{
			message("Download failed")
			return(NULL)
		}

		message("Formatting")
		o <- try(self$format_dataset())
		if('try-error' %in% class(o))
		{
			message("Formatting failed")
			message("Add ", self$ebi_id, " to ignore list")
			return(self$ebi_id)
		}

		message("Getting metadata")
		o <- try(self$organise_metadata())
		if('try-error' %in% class(o))
		{
			message("Formatting failed")
			return(NULL)
		}

		message("Upload metadata")
		o <- try(self$api_metadata_upload())
		if('try-error' %in% class(o))
		{
			message("GWAS upload failed")
			return(NULL)
		} else if(o$status_code != 200) {
			message("GWAS upload failed")
			print(httr::content(o))
			return(NULL)
		}

		message("Upload GWAS data")
		o <- try(self$api_gwas_upload())
		if('try-error' %in% class(o))
		{
			message("GWAS upload failed")
			self$upload_delete()
			return(NULL)
		} else if(o$status_code != 201) {
			message("GWAS upload failed")
			print(httr::content(o))
			self$upload_delete()
			return(NULL)
		}

		message("Completed successfully")
	}
))
