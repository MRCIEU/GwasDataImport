Dataset <- R6::R6Class("Dataset", list(
	filename = NULL,
	igd_id = NULL,
	wd = NULL,
	gwas_out = NULL,
	nsnp_read = NULL,
	nsnp = NULL,
	metadata = NULL,
	metadata_file = NULL,
	datainfo = NULL,
	datainfo_file = NULL,
	params = NULL,

	#' @description
	#' Initialise
	#' @param ebi_id e.g. GCST004426
	#' @param wd=tempdir() <what param does>
	#' @param ftp_path=NULL <what param does>
	#' @param igd_id Defaults to "ebi-a-<ebi_id>"
	#'
	#' @return new ObtainEbiDataset object
	initialize = function(filename, wd=tempdir(), igd_id=NULL)
	{
		if(!is.null(igd_id))
		{
			self$check_id(igd_id)
		}
		self$igd_id <- igd_id
		self$set_wd(wd)
		stopifnot(file.exists(filename))
		self$filename <- filename
	},

	check_id = function(id)
	{
		message("Checking existing IDs")
		r <- ieugwasr::api_query(paste0("gwasinfo/", id), access_token=ieugwasr::check_access_token(), method="GET")
		if(length(ieugwasr::get_query_content(r)) != 0)
		{
			stop("ID already in database: ", id)
		}
		message("Checking in-process IDs")
		r <- ieugwasr::api_query(paste0("edit/check/", id), access_token=ieugwasr::check_access_token(), method="GET")
		if(r$status_code == 404)
		{
			stop("Are you authenticated?")
		}
		if(length(ieugwasr::get_query_content(r)) != 0)
		{
			stop("ID already being processed: ", id)
		}
		return(NULL)
	},


	finalize = function()
	{
		message("Removing downloaded files")
		self$delete_wd()
	},

	#' @description
	#' delete working directory
	delete_wd = function()
	{
		message("Deleting download directory")
		unlink(self$wd, recursive=TRUE)
	},

	#' @description
	#' set working directory (creates)
	#' @param wd working directory
	set_wd = function(wd)
	{
		self$wd <- wd
		dir.create(self$wd, recursive=TRUE, showWarnings=FALSE)
	},

	determine_columns = function(params, nrows=100, gwas_file=self$filename, ...)
	{
		required_columns <- c("chr_col", "pos_col", "ea_col", "oa_col", "beta_col", "se_col", "pval_col")
		optional_columns <- c("snp_col", "eaf_col", "oaf_col", "ncase_col", "imp_z_col", "imp_info_col", "ncontrol_col")
		column_identifiers <- c(required_columns, optional_columns)

		stopifnot(all(required_columns %in% names(params)))

		m <- list(chr_col=1, pos_col=2, ea_col=3, oa_col=4, beta_col=5, se_col=6, pval_col=7)
		a <- data.table::fread(gwas_file, nrows=nrows, ...)
		if(is.infinite(nrows))
		{
			self$nsnp_read <- nrow(a)
		}
		out <- dplyr::tibble(
			chr=a[[params$chr_col]],
			pos=a[[params$pos_col]],
			ea=a[[params$ea_col]],
			oa=a[[params$oa_col]],
			beta=a[[params$beta_col]],
			se=a[[params$se_col]],
			pval=a[[params$pval_col]]
			) %>%
			subset(., ! (is.na(chr) | is.na(pos) | is.na(ea) | is.na(oa) | is.na(beta) | is.na(se) | is.na(pval)))
		stopifnot(all(is.numeric(out$pos)))
		stopifnot(all(is.numeric(out$beta)))
		stopifnot(all(is.numeric(out$se)))
		stopifnot(all(is.numeric(out$pval)))

		j <- length(required_columns)
		for(i in optional_columns)
		{
			if(!is.null(params[[i]]))
			{
				j <- j+1
				out[[gsub("_col$", "", i)]] <- a[[ params[[i]] ]]
				m[[i]] <- j
			}
		}
		stopifnot(nrow(out) > 0)
		m[["id"]] <- self$igd_id
		m[["header"]] <- "True"
		m[["gzipped"]] <- "True"
		m[["delimiter"]] <- "space"
		self$datainfo <- m
		self$params <- params
		message("Is this how the dataset should look:")
		print(str(out))
		if(is.infinite(nrows))
		{
			return(out)
		} else {
			return(NULL)
		}
	},


	#' @description
	#' format dataset
	#'
	#' @param dl=self$dl <what param does>
	format_dataset = function(gwas_file=self$filename, gwas_out = file.path(self$wd, "format.txt.gz"), params=self$params, ...)
	{
		out <- self$determine_columns(gwas_file=gwas_file, params=params, nrows=Inf, ...)
		self$datainfo[["gwas_file"]] <- gwas_out
		message("Determining build")
		out <- liftover_gwas(out)
		stopifnot(nrow(out) > 0)

		zz <- gzfile(gwas_out, "w")
		write.table(out, file=zz, row=F, col=TRUE, qu=FALSE, na="")
		close(zz)

		self$nsnp <- nrow(out)
		self$gwas_out <- gwas_out
		self$datainfo$gwas_file <- self$gwas_out
	},

	get_required_fields = function()
	{
		swagger <- jsonlite::read_json(paste0(options()$ieugwasr_api, "swagger.json"))
		p <- swagger$paths$"/edit/add"$post$parameters
		re <- sapply(p, function(x) "required" %in% names(x))
		re[re] <- sapply(p[re], function(x) x$required)
		dplyr::tibble(parameter=sapply(p, function(x) x$name), required=re) %>% 
		subset(., parameter != "X-Api-Token") %>%
		return
	},

	collect_metadata = function(metadata, igd_id=self$igd_id)
	{
		fields <- self$get_required_fields()
		stopifnot(all(fields$parameter[fields$required] %in% names(metadata)))
		l <- list()
		for(i in fields$parameter)
		{
			if(i %in% names(metadata))
			{
				l[[i]] <- metadata[[i]]
			}
		}
		l[["id"]] <- igd_id
		l[["nsnp_read"]] <- self$nsnp_read
		l[["nsnp"]] <- self$nsnp
		self$metadata <- l
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

	upload_metadata = function(metadata_file=self$metadata_file, access_token=ieugwasr::check_access_token())
	{
		ieugwasr::api_query("edit/add", query=self$metadata %>% magrittr::extract(keep_fields), access_token=access_token, method="POST")
	},

	upload_check = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		ieugwasr::api_query(paste0("edit/check/", id), access_token=access_token, method="GET")
	},

	upload_delete = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		ieugwasr::api_query(paste0("edit/delete/", id), access_token=access_token, method="DELETE")
	},

	upload_gwas = function(datainfo_file=self$datainfo_file, gwasfile=self$datainfo$filename, access_token=ieugwasr::check_access_token())
	{
		y <- x$datainfo
		y$gwas_file <- httr::upload_file(x$gwas_out)
		ieugwasr::api_query("edit/upload", query=y, access_token=access_token, method="POST", encode="multipart", timeout=600)
	}

))

