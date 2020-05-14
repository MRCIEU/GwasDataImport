ObtainEBIDataset <- R6::R6Class("ObtainEBIDataset", list(

	ebi_id = NULL,
	wd = NULL,
	dl = NULL,
	gwas_out = NULL,
	ftp_path = NULL,
	nsnp_read = NULL,
	nsnp = NULL,
	or_flag = NULL,
	metadata = NULL,
	metadata_file = NULL,
	datainfo = NULL,
	datainfo_file = NULL,

	initialize = function(ebi_id, wd=tempdir(), ftp_path=NULL)
	{
		self$ebi_id <- ebi_id
		self$set_wd(wd)
		self$ftp_path <- ftp_path
	},

	delete_wd = function()
	{
		message("Deleting download directory")
		unlink(self$wd, recursive=TRUE)
	},

	set_wd = function(wd)
	{
		self$wd <- wd
		dir.create(self$wd, recursive=TRUE)
	},

	download_dataset = function(ftp_path=self$ftp_path, ftp_url=options()$ebi_ftp_url, outdir=self$wd)
	{
		b <- basename(ftp_path)
		dl <- file.path(outdir, b)
		ftp <- file.path(ftp_url, ftp_path)
		cmd <- glue::glue("wget -O {dl} {ftp}")
		system(cmd)
		self$dl <- dl
	},

	format_dataset = function(dl=self$dl)
	{
		keep_cols <- c("hm_rsid", "hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value")

		a <- data.table::fread(dl, header=TRUE)
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
			dplyr::select(keep_cols) %>%
			subset(., !is.na(hm_chrom) & !is.na(hm_pos) & !is.na(hm_effect_allele) & !is.na(hm_beta) & !is.na(standard_error))

		stopifnot(nrow(out) > 0)

		gwas_out <- paste0(dl, ".format.gz")
		zz <- gzfile(gwas_out, "w")
		write.table(out, file=zz, row=F, col=TRUE, qu=FALSE, na="")
		close(zz)

		self$nsnp_read <- nrow(a)
		self$nsnp <- nrow(out)
		self$gwas_out <- gwas_out
		self$or_flag <- or_flag
	},

	organise_metadata = function(ebi_id=self$ebi_id, or_flag=self$or_flag, igd_id=NULL, units=NULL, sex="NA", category="NA", subcategory="NA", build="HG19/GRCh37", group_name="public")
	{
		l <- list()
		j <- jsonlite::read_json(paste0(options()$ebi_api, ebi_id))

		if(is.null(igd_id))
		{
			l[["id"]] <- paste0("ebi-a-", j[["accessionId"]])
		} else {
			l[["id"]] <- igd_id
		}
		l[["trait"]] <- j[["diseaseTrait"]][["trait"]]
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

	write_metadata = function(metadata=self$metadata, datainfo=self$datainfo, outdir=self$wd)
	{
		outfile <- file.path(outdir, "metadata.json")
		jsonlite::write_json(metadata, outfile, auto_unbox=TRUE)
		self$metadata_file <- outfile

		outfile <- file.path(outdir, "datainfo.json")
		jsonlite::write_json(datainfo, outfile, auto_unbox=TRUE)
		self$datainfo_file <- outfile
	},

	upload_metadata = function(metadata_file=self$metadata_file, datainfo_file=self$datainfo_file, access_token=ieugwasr::check_access_token())
	{
		keep_fields <- c("id", "trait", "note", "pmid", "year", "author", "population", "sample_size", "ncase", "ncontrol", "unit", "build", "group_name", "nsnp", "category", "subcategory", "sex")
		headers <- httr::add_headers(`X-Api-Token` = access_token, `X-Api-Source` = "EbiDataImport")
		httr::POST(
			paste0(options()$igd_api, "edit/add"),
			body = self$metadata %>% magrittr::extract(keep_fields),
			headers,
			encode = "json",
			httr::timeout(300)
		)
	},

	upload_check = function(id=self$metadata$id, access_token=ieugwasr::check_access_token())
	{
		headers <- httr::add_headers(`X-Api-Token` = access_token, `X-Api-Source` = "EbiDataImport")
		# httr::GET(
		# 	paste0(options()$igd_api, "edit/check/", id),
		# 	headers,
		# 	httr::timeout(300)
		# ) %>% httr::content()
		ieugwasr::editcheck(id)
	},

	upload_delete = function(id=self$metadata$id, access_token=ieugwasr::check_access_token())
	{
		headers <- httr::add_headers(`X-Api-Token` = access_token, `X-Api-Source` = "EbiDataImport")
		httr::DELETE(
			paste0(options()$igd_api, "edit/delete/", id),
			headers,
			httr::timeout(300)
		)
	},

	upload_gwas = function(datainfo_file=self$datainfo_file, gwasfile=self$datainfo$filename, access_token=ieugwasr::check_access_token())
	{
		headers <- httr::add_headers(`X-Api-Token` = access_token, `X-Api-Source` = "EbiDataImport")
		y <- x$datainfo
		y$gwas_file <- httr::upload_file(x$gwas_out)
		httr::POST(
			paste0(options()$igd_api, "edit/upload"),
			headers,
			body = y,
			httr::timeout(300)
		)
	}

))
