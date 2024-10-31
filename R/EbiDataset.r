#' Object that downloads, develops and uploads EBI dataset
#'
#' @export
EbiDataset <- R6::R6Class("EbiDataset", inherit = Dataset, list(

	#' @field ebi_id EBI ID to look for
	ebi_id = NULL,

	#' @field traitname Name of trait
	traitname = NULL,

	#' @field ftp_path Path to files in EBI FTP
	ftp_path = NULL,

	#' @field or_flag TRUE/FALSE if had to convert OR to beta
	or_flag = NULL,

	#' @field gwas_out1 Path to first look at EBI dataset
	gwas_out1 = NULL,

	#' @description
	#' Initialise object
	#' @param ebi_id e.g. GCST005522
	#' @param wd Directory in which to download and develop dataset. Default=tempdir(). Deleted automatically upon object removal
	#' @param ftp_path Pre-specified path to data. Default=NULL
	#' @param igd_id Defaults to "ebi-a-<ebi_id>"
	#' @param traitname Option to provide traitname of dataset
	#'
	#' @return A new EbiDataset object
	initialize = function(ebi_id, wd=tempdir(), ftp_path=NULL, igd_id=paste0("ebi-a-", ebi_id), traitname=NULL)
	{
		self$ebi_id <- ebi_id
		self$igd_id <- igd_id
		self$is_new_id(igd_id)
		self$set_wd(wd)
		self$ftp_path <- ftp_path
		self$traitname <- traitname
	},


	#' @description
	#' Download
	#' @param ftp_path  Pre-specified path to data. Default=self$ftp_path
	#' @param ftp_url Default=options()$ebi_ftp_url
	#' @param outdir Default=self$wd
	#' @importFrom glue glue
	download_dataset = function(ftp_path=self$ftp_path, outdir=self$wd, dl=TRUE)
	{
		dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
		b <- basename(ftp_path)
		filename <- file.path(outdir, b)
		ftp <- file.path(ftp_path)
		cmd <- glue::glue("wget -q -O {filename} {ftp}")
		if(dl) {
			system(cmd)
		}
		self$filename <- filename
	},


	#' @description
	#' organise data before formatting. This is slow but doesn't really matter
	#' @param filename Filename of GWAS dataset
	#' @param output Where to save formatted dataset
	format_ebi_dataset = function(filename=self$filename, output=file.path(self$wd, "step1.txt.gz"))
	{
		keep_cols <- c("hm_chrom", "hm_rsid", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value")

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
			dplyr::select(keep_cols)
		zz <- gzfile(output, "w")
		write.table(out, file=zz, row=F, col=TRUE, qu=FALSE)
		close(zz)
		self$gwas_out1 <- output
		self$or_flag <- or_flag
	},

	#' @description
	#' Download and parse metadata
	#' @param ebi_id Default=self$ebi_id
	#' @param or_flag Default=self$or_flag
	#' @param igd_id Default=NULL
	#' @param units Default=NULL
	#' @param sex Default="NA"
	#' @param category Default="NA"
	#' @param subcategory Default="NA"
	#' @param build Default="HG19/GRCh37"
	#' @param group_name Default="public"
	#' @param traitname Default=self$traitname
	organise_metadata = function(ebi_id=self$ebi_id, or_flag=self$or_flag, igd_id=self$igd_id, units=NULL, sex="NA", category="NA", subcategory="NA", build="HG19/GRCh37", group_name="public", traitname=self$traitname)
	{
		l <- list()
		# j <- jsonlite::read_json(paste0(options()$ebi_api, ebi_id))
		j <- gwasrapidd::get_studies(ebi_id)

		l[["id"]] <- igd_id
		if(is.null(traitname))
		{
			l[["trait"]] <- j@studies[["reported_trait"]]
		} else {
			l[["trait"]] <- traitname
		}
		l[["note"]] <- ""
		if(or_flag) l[["note"]] <- paste0(l[["note"]], "beta+se converted from OR+CI; ")
		if(!is.na(j@studies[["study_design_comment"]])) {
			l[["note"]] <- paste0(l[["note"]], j@studies[["study_design_comment"]], "; ")
		}
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
		l[["ontology"]] <- httr::GET(j[["_links"]][["efoTraits"]][["href"]]) %>% httr::content(., encoding="text") %>% 
			{.[["_embedded"]]} %>% 
			{.[["efoTraits"]]} %>%
			sapply(., function(x) x[["shortForm"]]) %>% paste(., collapse=";")
		l[["sex"]] <- sex
		l[["mr"]] <- 1
		l[["category"]] <- category
		l[["subcategory"]] <- subcategory
		self$metadata <- l

		self$datainfo[["cohort_cases"]] <- l[["ncase"]]
		self$datainfo[["cohort_controls"]] <- ifelse(is.na(l[["ncase"]]), l[["sample_size"]], l[["ncontrol"]])
	},

	#' @description
	#' Once initialised this function will string together everything i.e. downloading, processing and uploading
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
		o <- try({
			self$format_ebi_dataset()
			self$determine_columns(params=list(chr_col=1, snp_col=2, pos_col=3, oa_col=4, ea_col=5, eaf_col=6, beta_col=7, se_col=8, pval_col=9), gwas_file=x$gwas_out1)
			self$format_dataset(gwas_file=x$gwas_out1)
		})
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
		o <- try(self$api_gwasdata_upload())
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
