---
title: "Importing datasets into OpenGWAS [NEW WORKFLOW]"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Importing datasets into OpenGWAS [NEW WORKFLOW]}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Importing datasets into OpenGWAS [NEW WORKFLOW]

[[_TOC_]]

## As a Contributor

This role enables you to contribute datasets (create metadata, upload file for QC, check QC report and submit dataset for approval) to OpenGWAS.

You may contact Admin (details TBC) and request to be added as a Contributor. You will be granted access to https://api.opengwas.io/contribution.

You also need R/GwasDataImport installed and your OpenGWAS JWT (token) set up in your R environment.

## Setup

You can track your contributions via the web interface, but you will need R/GwasDataImport anyway to upload the files.

```{r}
library(GwasDataImport)

# Make sure you have switched to the "public" API base URL
ieugwasr::select_api()

# Make sure you are authenticated and have got the "contributor" or "admin" role 
ieugwasr::user()
```

## Metadata and dataset

For each dataset you will need to go through step 1 to 4. See also: [Import in bulk](#import-in-bulk).

### 1. Create metadata and obtain GWAS ID (two methods)

For each dataset to be uploaded, you can **either** use the web interface **or** the R/GwasDataImport package to create the metadata.

#### 1a. Use the web interface

https://api.opengwas.io/contribution/ provides a webform with dropdowns and tooltips for each field. This is handy when you are new to this process and/or only have a few datasets to upload.

GWAS ID will be available on the web interface for each dataset. Assume it's `ieu-b-5137`.

In R, specify the path (assume it's `~/bmi_test.txt.gz`) and GWAS ID, and then mark the metadata as uploaded.

```{r}
x <- Dataset$new(filename="~/bmi_test.txt.gz", igd_id='ieu-b-5137')

x$metadata_uploaded <- TRUE
```

#### 1b. Use R/GwasDataImport

**Alternatively**, you can opt to use R if you have an array of candidate datasets since it's easier to upload metadata of multiple datasets in a programmable way.
 
Assume the full path to the file is `~/bmi_test.txt.gz`.

```{r}
x <- Dataset$new(filename="~/bmi_test.txt.gz")

x$collect_metadata(list(
    trait="TEST - DO NOT USE 2",
    group_name="public",
    build="HG19/GRCh37",
    category="Risk factor",
    subcategory="Anthropometric",
    ontology="NA",
    population="Mixed",
    sex="Males and Females",
    sample_size=339224,
    author="Mendel GJ",
    year=2022,
    unit="SD"
))

x$api_metadata_upload()
```

GWAS ID will be returned by the last command. Assume it's `ieu-b-5137`. At the same time a new record will show up on https://api.opengwas.io/contribution/.

### 2. Modify the metadata (only when necessary)

You can modify the metadata either via the web interface (recommended) or through R/GwasDataImport, regardless of how you created the metadata (i.e. metadata created via the R package can be modified on the web interface, or vice versa).

Note that the metadata can only be modified when there is no QC pipeline associated, because at the last step when generating the QC report, the metadata will be hardcoded in the report. Metadata can only be modified when

- (a) the file has not been uploaded for QC, or
- (b) the QC pipeline had already finished and you decided to "delete uploaded files and QC products but keep the metadata", which effectly reverted the state to (a)

### 3. Format the file and upload for QC (R/GwasDataImport)

Always check that the GWAS ID and path information stored are accurate.

```{r}
x$igd_id
x$filename
```

Specify column mapping (1-indexed):

```{r}
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
```

Format the dataset (this may take a while):

```{r}
x$format_dataset()
```

Use the output to double-check the mapping. If necessary, run `x$determine_columns(...)` again with the correct mapping.

Upload the dataset (this may take a while):

```{r}
x$api_gwasdata_upload()
```

You will see the "Dataset has been added to the pipeline" message if the upload was successful.

And finally don't forget to clean up:

```{r}
x$delete_wd()
```

### 4. Check QC pipeline state and report, and submit for approval

On https://api.opengwas.io/contribution you can click the **2. QC** tab of the dataset popup and check pipeline state.

For each dataset, you should review the QC report when it's available and decide whether to submit the dataset for approval or not. You will have the following options:

- Submit the dataset for approval
    - Go to the **3. Approval & Release** tab and submit
- Re-upload the file for QC
    - Go to the **1. Metadata** tab and "delete uploaded files and QC products"
- Discard the dataset
    - Go to the **1. Metadata** tab and "delete everything"

You may also use the checkboxes on the main screen to select datasets and submit for approval in bulk.

## Import in bulk

If you have multiple datasets you may want to write an R snippet to semi-automate this process.

[Set up](#setup) for only once, then for each dataset go through [1b](#1b-use-rgwasdataimport), [2](#2-modify-the-metadata-only-when-necessary) and [3](#3-format-the-file-and-upload-for-qc-rgwasdataimport). Finally visit the portal in [4](#4-check-qc-pipeline-state-and-report-and-submit-for-approval) and use the checkboxes to submit in bulk.

## What's next

Admins will review and approve/reject each datasets.
