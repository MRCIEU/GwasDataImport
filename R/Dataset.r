#' Object that downloads, develops and uploads GWAS summary datasets for IEU OpenGWAS database
#'
#' @export
Dataset <- R6::R6Class("Dataset", list(
	
	#' @field filename Path to raw GWAS summary dataset
	filename = NULL,
	#' @field igd_id ID to use for upload
	igd_id = NULL,
	#' @field wd Work directory in which to save processed files. Will be deleted upon completion
	wd = NULL,
	#' @field gwas_out path to processed summary file
	gwas_out = NULL,
	#' @field nsnp_read Number of SNPs read initially
	nsnp_read = NULL,
	#' @field nsnp Number of SNPs retained after reading
	nsnp = NULL,
	#' @field metadata List of meta-data entries
	metadata = NULL,
	#' @field metadata_file Path to meta-data json file
	metadata_file = NULL,
	#' @field datainfo List of GWAS file parameters
	datainfo = NULL,
	#' @field datainfo_file Path to datainfo json file
	datainfo_file = NULL,
	#' @field params Initial column identifiers specified for raw dataset
	params = NULL,
	#' @field metadata_uploaded TRUE/FALSE of whether the metadata has been uploaded
	metadata_uploaded = FALSE,
	#' @field gwasdata_uploaded TRUE/FALSE of whether the gwas data has been uploaded
	gwasdata_uploaded = FALSE,

	#' @description
	#' Initialise
	#' @param filename Path to raw GWAS summary data file
	#' @param wd Path to directory to use as the working directory. Will be deleted upon completion - best to keep as the default randomly generated temporary directory
	#' @param igd_id Option to provide a specified ID for upload. If none provided then will use the next ieu-a batch ID
	#'
	#' @return new ObtainEbiDataset object
	initialize = function(filename=NULL, wd=tempdir(), igd_id=NULL)
	{
		if(!is.null(igd_id))
		{
			self$is_new_id(igd_id)
		}
		self$igd_id <- igd_id
		self$set_wd(wd)
		self$filename <- filename
	},

	#' @description
	#' Check if the specified ID is unique within the database. It checks published GWASs and those currently being processed
	#' @param id ID to check
	is_new_id = function(id=self$igd_id)
	{
		stopifnot(!is.null(id))
		message("Checking existing IDs")
		r <- ieugwasr::api_query(paste0("gwasinfo/", id), access_token=ieugwasr::check_access_token(), method="GET")
		if(length(ieugwasr::get_query_content(r)) != 0)
		{
			message("ID already in database: ", id)
			invisible(TRUE)
		}
		message("Checking in-process IDs")
		r <- ieugwasr::api_query(paste0("edit/check/", id), access_token=ieugwasr::check_access_token(), method="GET")
		if(r$status_code == 404)
		{
			stop("Are you authenticated?")
		}
		if(length(ieugwasr::get_query_content(r)) != 0)
		{
			message("ID already being processed: ", id)
			invisible(TRUE)
		}
		invisible(FALSE)
	},

	#' @description
	#' Delete working directory
	delete_wd = function()
	{
		message("Deleting temporary directory")
		unlink(self$wd, recursive=TRUE)
	},

	#' @description
	#' Set working directory (creates)
	#' @param wd working directory
	set_wd = function(wd)
	{
		self$wd <- wd
		message("Creating direcotry at:")
		message(wd)
		dir.create(self$wd, recursive=TRUE, showWarnings=FALSE)
	},

	se_from_bp = function(beta, pval, minp = 1e-300)
	{
		pval <- pmax(pval, minp)
		z <- qnorm(pval/2, low=FALSE)
		return(abs(beta / z))
	},

	#' @description
	#' Specify which columns in the dataset correspond to which fields. 
	#' @param params List of column identifiers. Identifiers can be numeric position or column header name. Required columns are: c("chr_col", "pos_col", "ea_col", "oa_col", "beta_col", "se_col", "pval_col"). Optional columns are: c("snp_col", "eaf_col", "oaf_col", "ncase_col", "imp_z_col", "imp_info_col", "ncontrol_col"). 
	#' @param nrows How many rows to read to check that parameters have been specified correctly
	#' @param gwas_file Filename to read
	#' @param ... Further arguments to pass to data.table::fread in order to correctly read the dataset
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
			)
		stopifnot(all(is.numeric(out$pos)))
		stopifnot(all(is.numeric(out$beta)))
		stopifnot(all(is.numeric(out$pval)))
		index <- is.na(out$se)
		if(any(index))
		{
			message("Deriving se for ", sum(index), " out of ", length(index), " entries")
			out$se[index] <- self$se_from_bp(out$beta, out$pval)
		}
		is_all_1 <- isTRUE(all.equal(out$se, rep(1, length(out$se)), tol=0.05))
		stopifnot(!is_all_1)


		j <- length(required_columns)
		for(i in optional_columns)
		{
			if(!is.null(params[[i]]))
			{
				if(any(!is.na( a[[ params[[i]] ]] )))
				{
					j <- j+1
					out[[gsub("_col$", "", i)]] <- a[[ params[[i]] ]]
					m[[i]] <- j
				}
			}
		}
		out <- out %>%
			subset(., ! (is.na(chr) | is.na(pos) | is.na(ea) | is.na(oa) | is.na(beta) | is.na(se) | is.na(pval)))
		stopifnot(nrow(out) > 0)

		m[["id"]] <- self$igd_id
		m[["header"]] <- "True"
		m[["gzipped"]] <- "True"
		m[["delimiter"]] <- "space"
		self$datainfo <- m
		self$params <- params
		message("Is this how the dataset should look?")
		print(str(out))
		if(is.infinite(nrows))
		{
			invisible(out)
		} else {
			invisible(NULL)
		}
	},


	#' @description
	#' Process dataset ready for uploading. Determins build and lifts over to hg19/b37 if necessary.
	#'
	#' @param gwas_file GWAS filename
	#' @param gwas_out Filename to save processed dataset to
	#' @param params Column specifications (see determine_columns for more info)
	#' @param ... Further arguments to pass to data.table::fread in order to correctly read the dataset
	format_dataset = function(gwas_file=self$filename, gwas_out = file.path(self$wd, "format.txt.gz"), params=self$params, ...)
	{
		out <- self$determine_columns(gwas_file=gwas_file, params=params, nrows=Inf, ...)
		self$datainfo[["gwas_file"]] <- gwas_out
		message("Determining build")
		out <- liftover_gwas(out, chr_col=self$datainfo$chr_col, pos_col=self$datainfo$pos_col, snp_col=self$datainfo$snp_col, ea_col=self$datainfo$ea_col, oa_col=self$datainfo$oa_col)
		stopifnot(nrow(out) > 0)

		zz <- gzfile(gwas_out, "w")
		write.table(out, file=zz, row=F, col=TRUE, qu=FALSE, na="")
		close(zz)

		message("Calculating md5")
		self$datainfo$md5 <- tools::md5sum(gwas_out)
		self$nsnp <- nrow(out)
		self$gwas_out <- gwas_out
		self$datainfo$gwas_file <- self$gwas_out
	},


	#' @description
	#' View the specifications for available meta data fields, as taken from http://gwas-api.mrcieu.ac.uk/docs 
	view_metadata_options = function()
	{
		swagger <- jsonlite::read_json(paste0(options()$ieugwasr_api, "swagger.json"))
		p <- swagger$paths$"/edit/add"$post$parameters
		str(p)
	},

	#' @description
	#' Get a list of fields and whether or not they are required
	#' @return data.frame
	get_required_fields = function()
	{
		swagger <- jsonlite::read_json(paste0(options()$ieugwasr_api, "swagger.json"))
		p <- swagger$paths$"/edit/add"$post$parameters
		re <- sapply(p, function(x) "required" %in% names(x))
		re[re] <- sapply(p[re], function(x) x$required)
		dplyr::tibble(parameter=sapply(p, function(x) x$name), required=re) %>% 
		subset(., parameter != "X-Api-Token") %>% 
		return()
	},

	#' @description
	#' Input metadata
	#' @param metadata List of meta-data fields and their values, see view_metadata_options for which fields need to be inputted.
	#' @param igd_id ID to be used for uploading to the database
	collect_metadata = function(metadata, igd_id=self$igd_id)
	{
		fields <- x$get_required_fields()
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
		excl <- names(metadata)[!names(metadata) %in% names(l)]
		if(length(excl) > 0)
		{
			message("The following field names were not permitted:")
			message(paste(excl, collapse="\n"))
		}
		self$metadata <- l
	},

	#' @description
	#' Write meta data to json file
	#' @param metadata List of meta data fields and their values
	#' @param datainfo List of data column parameters
	#' @param outdir Output directory to write json files
	write_metadata = function(metadata=self$metadata, datainfo=self$datainfo, outdir=self$wd)
	{
		outfile <- file.path(outdir, "metadata.json")
		jsonlite::write_json(metadata, outfile, auto_unbox=TRUE)
		self$metadata_file <- outfile

		outfile <- file.path(outdir, "datainfo.json")
		jsonlite::write_json(datainfo, outfile, auto_unbox=TRUE)
		self$datainfo_file <- outfile
	},

	#' @description
	#' Upload meta data to API
	#' @param metadata List of meta data fields and their values
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_metadata_upload = function(metadata=self$metadata, access_token=ieugwasr::check_access_token())
	{
		o <- ieugwasr::api_query("edit/add", query=metadata, access_token=access_token, method="POST")
		if(httr::status_code(o) == 200)
		{
			self$igd_id <- httr::content(o)$id
			self$metadata_uploaded <- TRUE
			message("ID: ", self$igd_id)
			message("Successfully uploaded metadata")
		} else {
			message("Failed to uploaded metadata")
		}
		return(o)
	},


	#' @description
	#' Upload meta data to API
	#' @param metadata List of meta data fields and their values
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_metadata_edit = function(metadata=self$metadata, access_token=ieugwasr::check_access_token())
	{
		stopifnot(!is.null(metadata$id))
		o <- ieugwasr::api_query("edit/edit", query=metadata, access_token=access_token, method="POST")
		if(httr::status_code(o) == 200)
		{
			self$igd_id <- httr::content(o)$id
			self$metadata_uploaded <- TRUE
			message("ID: ", self$igd_id)
			message("Successfully uploaded metadata")
		} else {
			message("Failed to uploaded metadata")
		}
		return(o)
	},

	#' @description
	#' View meta-data
	#' @param id ID to check
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_metadata_check = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		ieugwasr::api_query(paste0("edit/check/", id), access_token=access_token, method="GET")
	},

	#' @description
	#' Delete a dataset. This deletes the metadata AND any uploaded GWAS data (and related processing files)
	#' @param id ID to delete
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_metadata_delete = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		o <- ieugwasr::api_query(paste0("edit/delete/", id), access_token=access_token, method="DELETE")
		if(httr::status_code(o) == 200)
		{
			message("Successfully deleted gwas data and metadata from API")
			self$gwasdata_uploaded <- FALSE
			self$metadata_uploaded <- FALSE
		} else {
			message("Failed to delete gwas / meta data")
		}
	},

	#' @description
	#' Upload gwas dataset
	#' @param datainfo List of data column parameters
	#' @param gwasfile Path to processed gwasfile
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_gwasdata_upload = function(datainfo=self$datainfo, gwasfile=self$gwas_out, access_token=ieugwasr::check_access_token())
	{
		stopifnot(!is.null(self$igd_id))
		stopifnot(self$metadata_uploaded)
		y <- datainfo
		y$id <- self$igd_id
		y$gwas_file <- httr::upload_file(gwasfile)
		o <- ieugwasr::api_query("edit/upload", query=y, access_token=access_token, method="POST", encode="multipart", timeout=600)
		if(httr::status_code(o) == 201)
		{
			message("Successfully uploaded GWAS data")
			self$gwasdata_uploaded <- TRUE
		} else {
			message("Failed to upload GWAS data")
		}
		return(o)
	},

	#' @description
	#' Check status of API processing pipeline
	#' @param id ID to check
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_gwasdata_check = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		ieugwasr::api_query(paste0("quality_control/check/", id), access_token=access_token)
	},

	#' @description
	#' Delete a dataset. This deletes the metadata AND any uploaded GWAS data (and related processing files)
	#' @param id ID to delete
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_gwasdata_delete = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		self$api_metadata_delete(id=id, access_token=access_token)
	},

	#' @description
	#' View the html report for a processed dataset
	#' @param id ID of report to view
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_report = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		o <- httr::content(x$api_gwasdata_check())
		if(!any(grepl("html", unlist(o))))
		{
			message("html report hasn't been generated yet")
			return(o)
		}
		url <- paste0(options()$ieugwasr_api, "quality_control/report/", id)
		browseURL(url)
	},

	#' @description
	#' Release a dataset 
	#' @param comments Optional comments to provide when uploading
	#' @param passed_qc True or False
	#' @param id ID to release
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_gwas_release = function(comments=NULL, passed_qc="True", id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		payload <- list(id=id, comments=comments, passed_qc=passed_qc)
		ieugwasr::api_query("quality_control/release", query=payload, access_token=access_token, method="POST")
	}
))

