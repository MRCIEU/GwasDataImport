listftp <- function(url, R=TRUE)
{
	out <- tempfile()
	if(R)
	{
		cmd <- paste0("ncftpls -1 -R ", url, " > ", out)
	} else {
		cmd <- paste0("ncftpls -1 ", url, " > ", out)
	}
	system(cmd)
	a <- scan(out, what=character())
	a <- gsub("^./", "", a)
	return(a)
}


# filelist <- listftp(url, FALSE)
# filelistr <- listftp(url, TRUE) %>% grep("harmonised", ., value=TRUE) %>% grep("h.tsv.gz", ., value=TRUE)

# library(dplyr)
# temp <- tibble(filelistr) %>% tidyr::separate(filelistr, sep="/", into=letters[1:3])
# table(duplicated(temp$a))
# subset(temp, duplicated(a))

# table(duplicated(temp$c))
# head(filelist)


# head(filelistr)


# lapply(alldir, function(x)
# {
# 	getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE)

# })


# out <- html(url)

# id <- "GCST006906"



# read_n_lines <- function(id, list, url, n)
# {
# 	grep(paste0(id, "/harmonised"), b, value=TRUE)
# }



#' Get harmonised file for specific EBI ID
#'
#' @param ebi_id EBI ID e.g. GCST000879
#'
#' @export
#' @return ftp path (excluding server)
get_ftp_path <- function(ebi_id)
{
	ftplist <- listftp(options()$ebi_ftp_url, R=FALSE) %>% grep(paste0(ebi_id, "$"), ., value=TRUE)
	hfile <- listftp(file.path(options()$ebi_ftp_url, ftplist)) %>% grep("harmonised", ., value=TRUE) %>% grep("h.tsv.gz", ., value=TRUE)
	return(hfile)
}
