#' listftp
#'
#' <full description>
#'
#' @param url <what param does>
#' @param recursive=TRUE <what param does>
#'
#' @export
#' @return
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
#' <full description>
#'
#' @param ebi_ftp_url=options()$ebi_ftp_url <what param does>
#' @param blacklist=NULL <what param does>
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
#' <full description>
#'
#' @param dat <what param does>
#'
#' @export
#' @return
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
#' <full description>
#'
#' @param ebi_ftp_url=options()$ebi_ftp_url <what param does>
#'
#' @export
#' @return
ebi_datasets <- function(ebi_ftp_url=options()$ebi_ftp_url)
{
	message("Listing EBI FTP, may take a few minutes")
	l <- listftp(ebi_ftp_url, recursive=TRUE)
	all_dats <- l %>% grep("/harmonised/", ., value=TRUE) %>% grep("h.tsv.gz$", ., value=TRUE)
	message("Formatting...")
	dat <- dplyr::tibble(path=all_dats) %>% 
		tidyr::separate(path, sep="/", into=c("p1", "p2", "p3"), remove=FALSE) %>% 
		dplyr::select(-p2)
	dat$ebi_id <- do.call(rbind, strsplit(dat$p1, split="_"))[,3]
	dat <- dplyr::select(dat, path, ebi_id)	
	return(dat)
}


#' Get harmonised file for specific EBI ID
#'
#' @param ebi_id EBI ID e.g. GCST000879
#' @param ebi_ftp_url=options()$ebi_ftp_url <what param does>
#'
#' @export
#' @return ftp path (excluding server)
get_ftp_path <- function(ebi_id, ebi_ftp_url=options()$ebi_ftp_url)
{
	ftplist <- listftp(ebi_ftp_url, recursive=FALSE) %>% grep(paste0(ebi_id, "$"), ., value=TRUE)
	hfile <- listftp(file.path(options()$ebi_ftp_url, ftplist)) %>% grep("harmonised", ., value=TRUE) %>% grep("h.tsv.gz", ., value=TRUE)
	return(hfile)
}
