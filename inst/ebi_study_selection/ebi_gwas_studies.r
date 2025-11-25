library(gwasrapidd)
library(dplyr)
library(ggplot2)
library(lubridate)

ebi_studies <- readRDS("ebi_studies_20251124.rds")


ebi_studies@studies %>% glimpse

table(ebi_studies@studies$gxe)
table(ebi_studies@studies$gxg)
table(ebi_studies@studies$qualifier)
table(ebi_studies@studies$pooled)
table(ebi_studies@studies$full_pvalue_set)

# What data to get

# pubmed_id
# publication_date
# study_id
# reported_trait
# initial_sample_size
# replicaation_sample_size
# gxe
# gxg
# snp_count
# imputed
# pooled
# study_design_comment
# full_pvalue_set

studies <- ebi_studies@studies
pubs <- ebi_studies@publications %>% 
  dplyr::select(study_id, pubmed_id, publication_date, title)
studies <- left_join(studies, pubs, by="study_id")

dim(studies)
head(studies)

table(table(studies$pubmed_id))
tab <- group_by(pubs, pubmed_id) %>%
  summarise(n=n(), title=first(title), year=year(ymd(first(publication_date))))

subset(tab, n <= 100 & n > 50) %>%
{paste(.$year, .$n, .$title)}

large_p_keep <- c(
"Comprehensive genomic analysis of dietary habits in UK Biobank identifies hundreds of genetic associations.",
"Complex genetic signatures in immune cells underlie autoimmunity and inform therapy.",
"Insights into the genetic architecture of the human face.",
"A cross-platform approach identifies genetic regulators of human metabolism and health.",
"Cerebrospinal fluid metabolomics identifies 19 brain-related phenotype associations.",
"Clinical laboratory test-wide association scan of polygenic scores identifies biomarkers of complex disease.",
"Shared heritability of human face and brain shape.",
"Common genetic associations between age-related diseases.",
"Large-scale GWAS of food liking reveals genetic determinants and genetic correlations with distinct neurophysiological traits.",
"Genome-wide study on 72,298 individuals in Korean biobank data for 76 traits.",
"Genotyping and population characteristics of the China Kadoorie Biobank.",
"Genetic architecture of the structural connectome.",
"Diversity and longitudinal records: Genetic architecture of disease associations and polygenic risk in the Taiwanese Han population.",
"New Genetic Loci Implicated in Cardiac Morphology and Function Using Three-Dimensional Population Phenotyping.",
"Genome scans of facial features in East Africans and cross-population comparisons reveal novel associations."
)

table(tab$n)

table(large_p_keep %in% tab$title)

tab_keep <- subset(tab, title %in% large_p_keep)
keep <- subset(tab, n <=50)
keep <- bind_rows(tab_keep, keep)
dim(keep)
sum(keep$n)

studies_keep <- subset(studies, pubmed_id %in% keep$pubmed_id)
str(studies_keep)


studies_keep %>% filter(
  imputed, full_pvalue_set, !gxe, !gxg
) %>% {summary(.$snp_count)}

plot(1:10)

table(studies_keep$qualifier)
table(studies_keep$pooled)




parse_sample_size <- function(x) {
  # Handle NA values
  if (is.na(x)) return(NULL)

  # Remove commas from numbers
  b <- gsub("(\\d),(?=\\d)", "\\1", x, perl = TRUE) %>%
    # Split sample components based on ','
    strsplit(", ") %>%
    unlist()

  # Keep components that contain the word 'cases'
  b1 <- grep("cases", b, value = TRUE)

  # Keep components that contain the word 'controls'
  b2 <- grep("controls", b, value = TRUE)

  # Keep components that have neither
  b3 <- b[!b %in% c(b1, b2)]

  parse_bit <- function(b, name) {
    if(length(b) == 0) return(NULL)
    samples <- b %>%
      sapply(function(y) {
        suppressWarnings({
          num <- strsplit(y, " ") %>%
              unlist() %>% 
              as.numeric() %>% 
              na.omit() %>% 
              first()
            if (length(num) == 0) return(0)
            return(num)
          })
        })
    tibble(n=samples, text=names(samples)) %>%
      rowwise() %>% 
      mutate(text=trimws(gsub(as.character(n), "", text))) %>%
      mutate(what=name)
  }

  out <- bind_rows(
    parse_bit(b1, "cases"),
    parse_bit(b2, "controls"),
    parse_bit(b3, "sample_size")
  )
  print(out)
  if (!is.data.frame(out)) {
    return(NULL)
  }
  out
}







parse_sample_size(ebi_studies@studies$initial_sample_size[1100])
o <- parse_sample_size(ebi_studies@studies$replication_sample_size[110000])
b <- "461 European ancestry cases"
ebi_studies@studies$initial_sample_size[1]


ebi_studies@studies$initial_sample_size[10000]
