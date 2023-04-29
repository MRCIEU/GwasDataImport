#' Object that downloads, develops and uploads GWAS summary datasets for IEU OpenGWAS database
#'
#' @export
#' @importFrom R6 R6Class
Dataset <- R6::R6Class("Dataset", list(
	
	#' @field filename Path to raw GWAS summary dataset
	filename = NULL,
	#' @field igd_id ID to use for upload. If NULL then the next available ID in batch ieu-b will be used automatically
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
	#' @field metadata_test List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats
	metadata_test = NULL,
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
	#' @field metadata_upload_status Response from server about upload process
	metadata_upload_status = NULL,
	#' @field gwasdata_upload_status Response from server about upload process
	gwasdata_upload_status = NULL,

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
			invisible(FALSE)
		}
		message("Checking in-process IDs")
		r <- ieugwasr::api_query(paste0("edit/check/", id), access_token=ieugwasr::check_access_token(), method="GET")
		if(r$status_code == 404)
		{
			stop("Are you authenticated?")
		}
		invisible(TRUE)
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

	#' @description
	#' Estimate standard error from beta and p-value
	#' @param beta Effect size
	#' @param pval p-value
	#' @param minp Minimum p-value cutoff default = 1e-300
	se_from_bp = function(beta, pval, minp = 1e-300)
	{
		pval <- pmax(pval, minp)
		z <- qnorm(pval/2, low=FALSE)
		return(abs(beta / z))
	},

	#' @description
	#' Specify which columns in the dataset correspond to which fields. 
	#' @param params List of column identifiers. Identifiers can be numeric position or column header name. Required columns are: c("chr_col", "pos_col", "ea_col", "oa_col", "beta_col", "se_col", "pval_col","rsid_col"). Optional columns are: c("snp_col", "eaf_col", "oaf_col", "ncase_col", "imp_z_col", "imp_info_col", "ncontrol_col"). 
	#' @param nrows How many rows to read to check that parameters have been specified correctly
	#' @param gwas_file Filename to read
	#' @param ... Further arguments to pass to data.table::fread in order to correctly read the dataset
	#' @importFrom data.table fread
	determine_columns = function(params, nrows=100, gwas_file=self$filename, ...)
	{
		required_columns <- c("chr_col", "pos_col", "ea_col", "oa_col", "beta_col", "se_col", "pval_col")
		optional_columns <- c("snp_col", "eaf_col", "oaf_col", "ncase_col", "imp_z_col", "imp_info_col", "ncontrol_col")
		column_identifiers <- c(required_columns, optional_columns)

		
		stopifnot(all(required_columns %in% names(params)))

		m <- list(chr_col=1, pos_col=2, ea_col=3, oa_col=4, beta_col=5, se_col=6, pval_col=7)
		a <- data.table::fread(gwas_file, nrows=nrows, ...)
		# a <- data.table::fread(gwas_file, nrows=nrows)

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
		# head(out)
		stopifnot(all(is.numeric(out$pos)))
		stopifnot(all(is.numeric(out$beta)))
		stopifnot(all(is.numeric(out$pval)))
		index <- is.na(out$se)
		if(any(index))
		{
			message("Deriving se for ", sum(index), " out of ", length(index), " entries")
			out$se[index] <- self$se_from_bp(out$beta[index], out$pval[index])
		}
		is_all_1 <- isTRUE(all.equal(out$se, rep(1, length(out$se)), tol=0.05))

		if(is_all_1)
		{
			stop("beta values appear to be z-scores")
		}

		index.beta<-!is.na(out$beta)
		if(all(out$beta[index.beta]>0)){
			warning("beta values appear to be odds ratios")
		}

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

		message("Checking alleles are in A/C/T/G/D/I")
		problems <- grepl("[^ACTGactgdiDI]", out$ea) | grepl("[^ACTGactgdiDI]", out$oa)
		message(sum(problems), " variants with disallowed characters")
		out <- out[!problems, ]
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
	#' @param metadata_test  List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats
	#' @param ... Further arguments to pass to data.table::fread in order to correctly read the dataset
	
	format_dataset = function(gwas_file=self$filename, gwas_out = file.path(self$wd, "format.txt.gz"), params=self$params,metadata_test=self$metadata_test, ...)
	{
		# gwas_file=x$filename
		# gwas_out = file.path(x$wd, "format.txt.gz")
		# params=x$params
		# metadata_test=x$metadata_test
		

		# Check that CheckSumStats tests passed
		a<-unlist(metadata_test$eaf_conflicts)
		if(a[names(a)=="test"] == "fail") stop("The CheckSumStats test of the effect allele frequency column failed")
		a<-unlist(metadata_test$gc_conflicts)
		if(a[names(a)=="test"] == "strong_fail") stop("The CheckSumStats of the effect allele column failed")
		if(a[names(a)=="test"] == "moderate_fail") warning("The CheckSumStats of the effect allele column partially failed")
		a<-unlist(metadata_test$false_positive_hits)
		if(a[names(a)=="test"] == "fail") warning("The CheckSumStats test of reported GWAS hits failed. The GWAS hits in the test dataset are not present in the GWAS catalog")		


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
	#' @importFrom jsonlite read_json
	view_metadata_options = function()
	{
		swagger <- jsonlite::read_json(paste0(options()$ieugwasr_api, "swagger.json"))
		p <- swagger$paths$"/edit/add"$post$parameters
		str(p)
	},

	#' @description
	#' Get a list of GWAS data fields and whether or not they are required
	#' @return data.frame
	get_gwasdata_fields = function()
	{
		swagger <- jsonlite::read_json(paste0(options()$ieugwasr_api, "swagger.json"))
		p <- swagger$paths$"/edit/upload"$post$parameters
		re <- sapply(p, function(x) "required" %in% names(x))
		re[re] <- sapply(p[re], function(x) x$required)
		dplyr::tibble(
			parameter=sapply(p, function(x) x$name), 
			required=re, 
			description=sapply(p, function(x) x$description)
		) %>% 
		subset(., parameter != "X-Api-Token") %>%
		return()
	},


	#' @description
	#' Get a list of metadata fields and whether or not they are required
	#' @return data.frame
	get_metadata_fields = function()
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

		fields <- self$get_metadata_fields()
		required_fields <- c(fields$parameter[fields$required], "build", "unit", "ontology", "sample_size", "author", "year") %>% unique()
		if(!all(required_fields %in% names(metadata)))
		{
			stop("The following required fields were not provided: \n", paste(required_fields[!required_fields %in% names(metadata)], collapse="\n"))
		}
		stopifnot(metadata$build == "HG19/GRCh37")
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
	#' Check that the reported effect allele and effect allele frequency columns are correct. 
	#' @param gwas_file Filename to read
	#' @param params column names from x$determine_columns(). Required columns are: c("snp_col", "ea_col", "oa_col", "eaf_col" ) 
	#' @param metadata metadata from x$collect_metadata()

	check_meta_data = function(gwas_file=self$filename,params=self$params,metadata=self$metadata)
	{
		library(CheckSumStats)
		out<-data.table::fread(self$filename,nrows=Inf)	
		# out<-data.table::fread(x$filename,nrows=Inf)	
		
		
		#CheckSumStats expects these column names:
		names(out)[names(out) == params$snp_col]<-"rsid"
		names(out)[names(out) == params$ea_col]<-"effect_allele"
		names(out)[names(out) == params$oa_col]<-"other_allele"
		names(out)[names(out) == params$eaf_col]<-"eaf"
		names(out)[names(out) == params$beta_col]<-"beta"
		names(out)[names(out) == params$se_col]<-"se"
		names(out)[names(out) == params$pval_col]<-"p"
		
		out$population<-self$metadata$population #checksumstats expects an ancestry column 
		# out$population<-x$metadata$population
	
		# allele frequency conflicts with 1000 genomes super populations implying incorrect effect allele frequency column
		af_tests<-NULL
		af_tests<-af_conflicts_function(out=out)


		self$metadata_test$eaf_conflicts <- af_tests[1]
		self$metadata_test$eaf_conflicts_plot <- af_tests[2]

		# x$metadata_test$eaf_conflicts <- af_tests[1]
		# x$metadata_test$eaf_conflicts_plot <- af_tests[2]

		
		# effect size conflicts with GWAS catalog implying incorrect effect allele column
		# this function can be slow if: 1) there are lots of "GWAS hits" for the trait of interest, 2) when searching on both EFO and trait name and 3) when mapping genetic associations to study_id in the GWAS catalog (this latter step is needed for matching of associations on ancestry). To speed things up, the default is to only search on trait name and to ignore study ancestry. If not associations are found on trait name, the search is then done on the user provided ontology. It is assumed that ontology is a required metadata field
		
		gc_tests<-NULL
		gc_tests<-gc_conflicts_function(out=out,metadata=metadata)	
		# efo_id=NULL,trait=trait,efo_id=NULL,ignore_conflict=FALSE,map_association_to_study=FALSE		
		self$metadata_test$gc_conflicts <- gc_tests[1]
		self$metadata_test$gc_conflicts_plot <- gc_tests[2]

		# find reported GWAS hits in the gWAS catalog. absence suggests problems with the provided summary data. e.g. failure to exclude false positives / unreliable associations / to pass GWAS results throigh post GWAS QC
		false_positive_tests<-NULL
		# x$metadata_test$eaf_conflicts
		false_positive_tests<-reported_gwas_hits_in_gwascatalog(out=out,eaf_test=af_tests[1],metadata=metadata)
		# false_positive_tests<-reported_gwas_hits_in_gwascatalog(out=out,eaf_test=x$metadata_test$eaf_conflicts,metadata=metadata)
		self$metadata_test$false_positive_hits<-false_positive_tests


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
	#' @param metadata_test  List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	#' @importFrom httr status_code
	api_metadata_upload = function(metadata=self$metadata, metadata_test =self$metadata_test, access_token=ieugwasr::check_access_token())
	{

		# check that CheckSumStats tests passed
		a<-unlist(self$metadata_test$eaf_conflicts)
		if(a[names(a)=="test"] == "fail") stop("effect allele frequency column incorrectly specified")
		a<-unlist(self$metadata_test$gc_conflicts)
		if(a[names(a)=="test"] == "strong_fail") stop("effect allele column incorrectly specified")
		if(a[names(a)=="test"] == "moderate_fail") warning("effect allele column may be incorrectly specified")
		a<-unlist(self$metadata_test$false_positive_hits)
		if(a[names(a)=="test"] == "fail") warning("The CheckSumStats test of reported GWAS hits failed. The GWAS hits in the test dataset are not present in the GWAS catalog")	

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
		self$metadata_upload_status <- o
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
		self$metadata_upload_status <- o
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
		self$metadata_upload_status <- o
	},

	#' @description
	#' Upload gwas dataset
	#' @param datainfo List of data column parameters
	#' @param gwasfile Path to processed gwasfile
	#' @param metadata_test  List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_gwasdata_upload = function(datainfo=self$datainfo, gwasfile=self$gwas_out,metadata_test=self$metadata_test, access_token=ieugwasr::check_access_token())
	{
		# check that CheckSumStats tests passed
		a<-unlist(self$metadata_test$eaf_conflicts)
		if(a[names(a)=="test"] == "fail") stop("effect allele frequency column incorrectly specified")
		a<-unlist(self$metadata_test$gc_conflicts)
		if(a[names(a)=="test"] == "strong_fail") stop("effect allele column incorrectly specified")
		if(a[names(a)=="test"] == "moderate_fail") warning("effect allele column may be incorrectly specified")
		a<-unlist(self$metadata_test$false_positive_hits)
		if(a[names(a)=="test"] == "fail") warning("The CheckSumStats test of reported GWAS hits failed. The GWAS hits in the test dataset are not present in the GWAS catalog")	

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
		self$gwasdata_upload_status <- o
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
	#' Check the status of the GWAS QC processing pipeline
	#' @param id ID to delete
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_qc_status = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		readr::read_file(paste0(
			options()$cromwell_api, "/api/workflows/v1/query?label=gwas_id:", self$igd_id)
			) %>% jsonlite::prettify()
	},

	#' @description
	#' View the html report for a processed dataset
	#' @param id ID of report to view
	#' @param access_token Google OAuth2.0 token. See ieugwasr documentation for more info
	api_report = function(id=self$igd_id, access_token=ieugwasr::check_access_token())
	{
		o <- httr::content(self$api_gwasdata_check())
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

af_conflicts_function<-function(out=NULL)
{
	af_conflicts<-CheckSumStats::flag_af_conflicts(target_dat=out)
	af_conflicts_test<-"test not possible"
	if(af_conflicts$number_of_snps>10 & af_conflicts$proportion_conflicts>=0.10){
		af_conflicts_test<-"fail"
	}		
	if(af_conflicts$number_of_snps>10 & af_conflicts$proportion_conflicts<0.10){
		af_conflicts_test<-"pass"
	}		
	# if(af_conflicts_test=="fail") warning("allele frequency conflicts identified. It looks like you have incorrectly specified the effect allele frequency column.")
	af_conflicts<-c(af_conflicts,"test"=af_conflicts_test)
	Plot1<-CheckSumStats::make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=out)
	# af_conflicts<-c(af_conflicts,"plot"=Plot1)
	return(list(af_conflicts,Plot1))
}

# stop("allele frequency conflicts identified. It looks like you have incorrectly specified the effect allele frequency column.")
# , efo_id=NULL,trait=NULL,ignore_conflict=FALSE,map_association_to_study=FALSE
gc_conflicts_function<-function(out=NULL,metadata=NULL)
{

	gwas_catalog<-CheckSumStats::gwas_catalog_hits(trait=metadata$trait)
	if(length(gwas_catalog) == 1)
	{
		if(gwas_catalog == "no results found")
		{
			gwas_catalog<-CheckSumStats::gwas_catalog_hits(efo_id=metadata$ontology)
		}
	}
	gc_dat<-CheckSumStats::compare_effect_to_gwascatalog2(dat=out,gwas_catalog=gwas_catalog,trait=metadata$trait)
	# if(gc_dat == "no results found"){
	# 	gc_dat<-CheckSumStats::compare_effect_to_gwascatalog2(dat=out,efo_id=metadata$ontology)
	# }

	gc_conflicts<-CheckSumStats::flag_gc_conflicts2(gc_dat=gc_dat)	
	
	N_snps<-gc_conflicts$effect_size_conflicts$n_snps
	
	n_high<-unlist(gc_conflicts$effect_size_conflicts["high conflict"])
	n_moderate<-unlist(gc_conflicts$effect_size_conflicts["moderate conflict"])
	N_conflicts<-n_high + n_moderate
	proportion_conflicts<-N_conflicts/N_snps

	gc_conflicts_test<-"Number of SNPs may be too low for reliable test. Interpret with caution"
	if(N_snps>10){
		if(proportion_conflicts >=0.8){ 
			gc_conflicts_test<-"strong_fail"
			# warning("Very strong effect size conflicts with GWAS catalog identified. The effect allele column looks wrong")
		}
		if(proportion_conflicts<0.8 & proportion_conflicts>0.3){
			gc_conflicts_test<-"moderate_fail"
			# warning("Effect size conflicts with GWAS catalog identified. It looks like you may have incorrectly specified the effect allele column")
		}
		if(proportion_conflicts<0.3){
			gc_conflicts_test<-"pass"
		}
	}
	gc_conflicts<-c(gc_conflicts,"test"=gc_conflicts_test)
	Plot2<-CheckSumStats::make_plot_gwas_catalog(dat=out,gwas_catalog=gwas_catalog,gc_dat=gc_dat)
	return(list(gc_conflicts,Plot2))
}

reported_gwas_hits_in_gwascatalog<-function(out=NULL,clump=TRUE,eaf_test=NULL,metadata=NULL)
{	
	Pop<-infer_ancestry(target_dat=out)
	l2<-unlist(Pop)
	Pop<-names(l2[which(l2==max(l2))])
	dat<-out[which(out$p<5e-8),]
	# l<-unlist(eaf_test)
	# if(l[names(l)=="test"] =="fail") warning("eaf test failed, which might undermine the infer_ancestry test for ld_clumping")
	if(nrow(dat)==0) return("no GWAS hits at 5e-8 threshold")

	Clump<-ieugwasr::ld_clump(dat=dplyr::tibble(rsid=dat$rsid, pval=dat$p),clump_r2 = 0.001,clump_p=5e-8,pop=Pop) #it seems like there is a bug such that when none of the variants are present in the reference panel this function returns a Server code:503.  
	gwas_hits<-Clump$rsid
	if(length(Clump$rsid)<5)
	{
		gwas_hits<-dat$rsid
	}

	gc_list<-find_hits_in_gwas_catalog(gwas_hits=gwas_hits,trait=metadata$trait,distance_threshold=50000) 
	#if no results are found on reported trait, then use the user provided ontology to find associations
	if(!is(gc_list, "list")){
		if(gc_list == "no results found")
		{
			metadata$ontology<-gsub(":","_",metadata$ontology) #currently there is a bug in the gwasrapidd code that does not allow EFOs with semicolons. The semicolon must be replaced with an underscore. I have writtent to the developers to see if this can be fixed.  
			gc_list<-find_hits_in_gwas_catalog(gwas_hits=gwas_hits,efo_id=metadata$ontology,distance_threshold=50000) 
		}
	}

	if(!is(gc_list, "list")){
		if(gc_list == "no results found") return(gc_list)
	}

	# metadata$ontology

	N_in_gc<-length(gc_list$in_gc)
	N_not_in_gc<-length(gc_list$not_in_gc )
	N<-N_not_in_gc+N_in_gc
	test_hits_gc<-"test not possible"
	if(N_not_in_gc/N>0.3)
	{
		test_hits_gc<-fail
	}
	if(N_not_in_gc/N<0.3)
	{
		test_hits_gc<-"pass"
	}
	false_positive_hits<-c(gc_list,"test"=test_hits_gc)

	return(false_positive_hits)
}
