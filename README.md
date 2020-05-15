
# EbiDataImport

<!-- badges: start -->
<!-- badges: end -->

The goal of EbiDataImport is to ...

## Installation

You can install the released version of EbiDataImport from [CRAN](https://CRAN.R-project.org) with:

``` r
remotes::install_github("mrcieu/EbiDataImport")
```


## Outline

1. Download and format the data and metadata for an EBI dataset
2. Upload the formatted EBI data and metadata to IGD API
3. Determine which datasets in EBI are not present in IGD. Allow a 'blacklist' of EBI datasets to ignore. (TODO)

See `vignettes` and `tests` for examples.








## Ignore list

Need to figure out a solution for storing this. Many EBI datasets don't have sufficient info to be included, so we need to ignore them instead of trying to upload them each time we look for new studies.

On the `ieu-ssd` machine where the original set of EBI datasets was downloaded and processed, this is how the first ignore list was determined:

```
library(dplyr)
a <- read.csv("/data/ebi_gwascatalog/pmids.txt")
b <- scan("/data/ebi_gwascatalog/study_ids.txt", character())
got <- ieugwasr::gwasinfo() %>% subset(., grepl("ebi-a", id)) %>% {.$id} %>% gsub("ebi-a-", "", .)
a$got <- a$StudyID %in% got
head(a)
table(a$got)
table(a$got, a$Status)
ebi_ignore_list <- a$StudyID[a$Status == "harmonised" & !a$got]
write.table(ebi_ignore_list, file="ebi_ignore_list.txt", row=F, col=F, qu=F)
```

