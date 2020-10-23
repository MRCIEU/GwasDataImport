.onAttach <- function(libname, pkgname) {
	options(ebi_ftp_url="ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/")
	options(ebi_api="https://www.ebi.ac.uk/gwas/rest/api/studies/")
	options(cromwell_api="http://ieu-db-interface.epi.bris.ac.uk:8001")
	# options(igd_api="http://ieu-db-interface.epi.bris.ac.uk:8082/")
	# options(igd_api="http://127.0.0.1:5000/")
	ieugwasr::select_api("private")
}
