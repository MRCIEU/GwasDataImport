#' listftp
#'
#' List all files on the EBI FTP server
#'
#' @param url FTP url to look up
#' @param recursive If false then just the top directory, otherwise list recursively
#'
#' @export
#' @return Vector of paths
listftp <- function(url=options()$ebi_ftp_url, recursive=TRUE)
{
	out <- tempfile()
	if(recursive)
	{
		cmd <- paste0("ncftpls -1 -R ", url, " > ", out)
	} else {
		cmd <- paste0("ncftpls -1 ", url, " > ", out)
	}
	system(cmd)
	a <- scan(out, what=character())
	a <- gsub("^./", "", a)
	unlink(out)
	return(a)
}


#' determine_new_datasets
#'
#' Figure out which datasets are not present in the database
#'
#' @param ebi_ftp_url FTP url default=options()$ebi_ftp_url
#' @param blacklist List of EBI datasets to ignore default=NULL
#' @param exclude_multi_datasets If a EBI ID has more than one dataset then should it be ignored
#'
#' @export
#' @return data frame
determine_new_datasets <- function(ebi_ftp_url=options()$ebi_ftp_url, blacklist=NULL, exclude_multi_datasets=TRUE)
{
	message("Getting IGD dataset list")
	igd_ids <- suppressMessages(ieugwasr::gwasinfo()) %>% subset(., grepl("^ebi-a", id)) %>% {gsub("^ebi-a-", "", .$id)}
	message("Found ", length(igd_ids), " datasets in IGD")

	message("Determining existing EBI datasets")
	dat <- ebi_datasets(ebi_ftp_url)

	if(exclude_multi_datasets)
	{
		dup_ids <- dat$ebi_id[duplicated(dat$ebi_id)] %>% unique
		dat <- subset(dat, !ebi_id %in% dup_ids)
	}

	dat$need <- ! dat$ebi_id %in% igd_ids
	if(!is.null(blacklist))
	{
		dat <- subset(dat, ! ebi_id %in% blacklist)
	}

	dat <- subset(dat, need)

	return(dat)
}


#' being_processed
#'
#' List of EBI datasets that are currently being processed
#'
#' @param dat Output from \code{determine_new_datasets}
#'
#' @export
#' @return Updated dat
being_processed <- function(dat)
{
	message("Checking if new datasets aren't in check list")
	for(i in 1:nrow(dat))
	{
		out <- suppressMessages(ieugwasr::editcheck(paste0("ebi-a-", dat$ebi_id[i])))
		if(length(out) != 0)
		{
			dat$need[i] <- FALSE
		}
		message(i, " of ", nrow(dat), ": ", dat$ebi_id[i], "; need: ", dat$need[i])
	}
	return(dat)	
}

#' ebi_datasets
#'
#' Convert output from listftp into something that is easier to read
#'
#' @param ebi_ftp_url EBI FTP default=options()$ebi_ftp_url
#'
#' @export
#' @return data frame
ebi_datasets <- function(ebi_ftp_url=options()$ebi_ftp_url)
{
	message("Listing EBI FTP, may take a few minutes")
	l <- listftp(ebi_ftp_url, recursive=TRUE)
	all_dats <- l %>% grep("/harmonised/", ., value=TRUE) %>% grep("h.tsv.gz$", ., value=TRUE)
	message("Formatting...")
	dat <- dplyr::tibble(path=all_dats) %>% 
		tidyr::separate(path, sep="/", into=c("p0", "p1", "p2", "p3"), remove=FALSE) %>% 
		dplyr::select(-c(p0, p2))
	dat$ebi_id <- dat$p1
	dat <- dplyr::select(dat, path, ebi_id)	
	return(dat)
}


#' Get harmonised file for specific EBI ID
#'
#' @param ebi_id EBI ID e.g. GCST000879
#' @param ebi_ftp_url EBI FTP default=options()$ebi_ftp_url
#'
#' @export
#' @return ftp path (excluding server)
get_ftp_path <- function(ebi_id, ebi_ftp_url=options()$ebi_ftp_url)
{
	ftplist <- listftp(ebi_ftp_url, recursive=FALSE) %>% grep(paste0(ebi_id, "$"), ., value=TRUE)
	hfile <- listftp(file.path(options()$ebi_ftp_url, ftplist)) %>% grep("harmonised", ., value=TRUE) %>% grep("h.tsv.gz", ., value=TRUE)
	return(hfile)
}
