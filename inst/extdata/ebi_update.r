# - Step 2 instructions - point them to the URL for the R package https://mrcieu.github.io/GwasDataImport/
# - R package to pull down existing meta data


get_ebi_list <- function() {
	harmonised <- scan(url(file.path(options()$ebi_ftp_url, "harmonised_list.txt")), what=character()) %>%
		gsub("^./", options()$ebi_ftp_url, .)
	
	head(harmonised)


	dat <- dplyr::tibble(
		path=harmonised, 
		fn=basename(path) %>% gsub(".h.tsv.gz", "", .),
		yaml=paste0(path, "-meta.yaml")
		) %>% 
		tidyr::separate(fn, sep="-", into=c("pubmedid", "ebi_id", "EFO"), remove=FALSE) 
	
	ind <- is.na(dat$ebi_id)
	dat$ebi_id[ind] <- dat$pubmedid[ind]
	dat$pubmedid[ind] <- NA

	library(yaml)
	library(purrr)

	yam <- purrr::map(dat$yaml, yaml::yaml.load_file, .progress=TRUE)

	yaml::yaml.load_file(dat$yaml[4])
	yam <- list()
	for(i in 2296:length(dat$yaml)) {
		message(i)
		yam[[i]] <- try(yaml.load_file(dat$yaml[i]))
	}

	dat$path[4]
	table(is.na(dat$ebi_id), is.na(dat$EFO))


	table(dat$ebi_id %in% got)
	table(got %in% dat$ebi_id)

	dim(dat)

	table(table(dat$pubmedid))

	table(!is.na(dat$pubmedid), dat$ebi_id %in% got)



	got <- ieugwasr::gwasinfo() %>% subset(., grepl("^ebi-a", id)) %>% {gsub("^ebi-a-", "", .$id)}
	datp <- subset(dat, !is.na(dat$pubmedid))
	table(table(datp$pubmedid))
	library(dplyr)
	group_by(datp, pubmedid) %>% summarise(n=n()) %>% arrange(-n) %>% subset(., n>100)

no - 4753 - 35078996 - decode proteomics
no - 3695 - 34662886 - backman exomes
no - 2939 - 34737426 - fastgwa
yes - 1863 - 35668104 - lipids in 5k samples
yes - 922 - 34503513 - metabolites in 5k south asians
yes - 731 - 32929287 - immune cell traits in 3.8k sardinians
yes - 469 - 35115689 - gut microbiome in 6k samples
no - 433 - 34594039 - cross pop gwas
yes - 430 - 33462482 - gut microbiome in 9k samples (German)
yes - 412 - 35115690 - gut microbiome in 8k samples (Dutch)
no - 337 - 33437055 - CSF metabolites in 291 samples
no - 257 - 33303764 - wgs proteomics in 1.4k samples
no - 249 - 35213538 - ukb metabolites 
yes - 211 - 33462485 - gut microbiome in 18k samples
no - 161 - 34017140 - regenie
yes - 134 - 33414548 - blood metabolites in 80k
no - 115 - 33959723 - disease gwas in ukb
no - 110 - 35264221 - plasma proteomics in 400 samples
no - 106 - 33328453 - plasma proteomics in 10k samples



datnop <- subset(dat, is.na(pubmedid))


library(gwasrapidd)

ebi_ids <- datnop$ebi_id
ebi_ids_s <- split(ebi_ids, 1:20)
str(ebi_ids_s)


ebi_ids_s[[1]][2328]
l <- list()
for(i in 2:length(ebi_ids_s)) {
	message(i)
	s <- ebi_ids_s[[i]]
	s <- grep("_", s, value=TRUE, invert=TRUE)
	length(s)
	l[[i]] <- gwasrapidd::get_studies(s)
}
a <- gwasrapidd::get_studies(datnop$ebi_id[1:100])
a
split(1:10, 1:3)

datnop$path


slotNames(l[[1]])
l[[1]]@studies
temp <- l[[1]]@publications

table(table(temp$pubmed_id))


l1 <- lapply(l, \(x) x@publications) %>% bind_rows()
l1
l1 %>% group_by(pubmed_id) %>% summarise(n=dplyr::n()) %>% arrange(-n) %>% subset(., n>500) %>% as.data.frame %>% {sum(.$n)}



   pubmed_id    n
1   38811844 6828 - brain imaging
2   35870639 6763 - pqtls
3   36635386 4443 - plasma metabolites
4   33875891 3841 - brain imaging (got)
5   29875488 3283 - plasma proteome
6   35995766 3037 - metabolites (black)
7   37277652 2673 - metabolites
8   36168886 1299 - pqtls
9   38412862 1221 - protein ratios
10  35120996 1185 - blood metabolites
11  24816252  486 - blood metabolites (got)
12  37524825  446 - eicosanoids
13  38565889  341 - plasma proteomics
14  39278973  325 - metabolites ukb
15  34662886  291 - backman exomes
16  38580839  279 - brain imaging
17  33821002  274
18  36764567  249
19  36349687  244
20  38448586  233
21  38438384  206
22  35585065  186
23  34857772  184
24  36261456  149
25  32193382  139
26  26068415  128
27  27005778  123
28  33288918  122
29  38351177  112


x <- Dataset$new(filename="35078996-GCST90086510-EFO_0007937.h.tsv.gz")

gwasrapidd::get_studies("GCST90086510")@studies %>% str

x$determine_columns(list(
  chr_col=11,
  pos_col=12,
  ea_col=2,
  oa_col=3,
  beta_col=4,
  se_col=5,
  pval_col=6,
  snp_col=1,
  eaf_col=7,
  ncontrol_col=8
))



harmonised <- scan(url(file.path(options()$ebi_ftp_url, "harmonised_list.txt")), what=character()) %>%
    gsub("^./", options()$ebi_ftp_url, .)
got <- ieugwasr::gwasinfo() %>% subset(., grepl("^ebi-a", id)) %>% {gsub("^ebi-a-", "", .$id)}
datp <- subset(dat, !is.na(dat$pubmedid))
excl <- c(35078996, 34662886, 34737426, 34594039, 33437055, 33303764, 35213538, 34017140, 33959723, 35264221, 33328453)
pub_to_keep <- datp %>% subset(., ! ebi_id %in% got & ! pubmedid %in% excl)


library(devtools)
load_all()

ebi_id <- pub_to_keep$ebi_id[1]
dat$path[dat$ebi_id==ebi_id]
wd <- x$wd
x <- EbiDataset$new(
    ebi_id = ebi_id, 
    ftp_path = dat$path[dat$ebi_id==ebi_id],
    igd_id = paste0("ebi-a-", ebi_id),
    wd = wd
)
x$download_dataset(dl=FALSE)
x$format_ebi_dataset()
x$determine_columns(params=list(chr_col=1, snp_col=2, pos_col=3, oa_col=4, ea_col=5, eaf_col=6, beta_col=7, se_col=8, pval_col=9), gwas_file=x$gwas_out1)
x$format_dataset(gwas_file=x$gwas_out1)
x$organise_metadata()

ieugwasr::api_query(paste0("edit/check/", ebi_id), method="GET") %>% httr::content() %>% length