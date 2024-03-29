% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dataset.r
\name{Dataset}
\alias{Dataset}
\title{Object that downloads, develops and uploads GWAS summary datasets for IEU OpenGWAS database}
\description{
Object that downloads, develops and uploads GWAS summary datasets for IEU OpenGWAS database

Object that downloads, develops and uploads GWAS summary datasets for IEU OpenGWAS database
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{filename}}{Path to raw GWAS summary dataset}

\item{\code{igd_id}}{ID to use for upload. If NULL then the next available ID in batch ieu-b will be used automatically}

\item{\code{wd}}{Work directory in which to save processed files. Will be deleted upon completion}

\item{\code{gwas_out}}{path to processed summary file}

\item{\code{nsnp_read}}{Number of SNPs read initially}

\item{\code{nsnp}}{Number of SNPs retained after reading}

\item{\code{metadata}}{List of meta-data entries}

\item{\code{metadata_test}}{List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats}

\item{\code{metadata_file}}{Path to meta-data json file}

\item{\code{datainfo}}{List of GWAS file parameters}

\item{\code{datainfo_file}}{Path to datainfo json file}

\item{\code{params}}{Initial column identifiers specified for raw dataset}

\item{\code{metadata_uploaded}}{TRUE/FALSE of whether the metadata has been uploaded}

\item{\code{gwasdata_uploaded}}{TRUE/FALSE of whether the gwas data has been uploaded}

\item{\code{metadata_upload_status}}{Response from server about upload process}

\item{\code{gwasdata_upload_status}}{Response from server about upload process}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Dataset-new}{\code{Dataset$new()}}
\item \href{#method-Dataset-is_new_id}{\code{Dataset$is_new_id()}}
\item \href{#method-Dataset-delete_wd}{\code{Dataset$delete_wd()}}
\item \href{#method-Dataset-set_wd}{\code{Dataset$set_wd()}}
\item \href{#method-Dataset-se_from_bp}{\code{Dataset$se_from_bp()}}
\item \href{#method-Dataset-determine_columns}{\code{Dataset$determine_columns()}}
\item \href{#method-Dataset-format_dataset}{\code{Dataset$format_dataset()}}
\item \href{#method-Dataset-view_metadata_options}{\code{Dataset$view_metadata_options()}}
\item \href{#method-Dataset-get_gwasdata_fields}{\code{Dataset$get_gwasdata_fields()}}
\item \href{#method-Dataset-get_metadata_fields}{\code{Dataset$get_metadata_fields()}}
\item \href{#method-Dataset-collect_metadata}{\code{Dataset$collect_metadata()}}
\item \href{#method-Dataset-check_meta_data}{\code{Dataset$check_meta_data()}}
\item \href{#method-Dataset-write_metadata}{\code{Dataset$write_metadata()}}
\item \href{#method-Dataset-api_metadata_upload}{\code{Dataset$api_metadata_upload()}}
\item \href{#method-Dataset-api_metadata_edit}{\code{Dataset$api_metadata_edit()}}
\item \href{#method-Dataset-api_metadata_check}{\code{Dataset$api_metadata_check()}}
\item \href{#method-Dataset-api_metadata_delete}{\code{Dataset$api_metadata_delete()}}
\item \href{#method-Dataset-api_gwasdata_upload}{\code{Dataset$api_gwasdata_upload()}}
\item \href{#method-Dataset-api_gwasdata_check}{\code{Dataset$api_gwasdata_check()}}
\item \href{#method-Dataset-api_gwasdata_delete}{\code{Dataset$api_gwasdata_delete()}}
\item \href{#method-Dataset-api_qc_status}{\code{Dataset$api_qc_status()}}
\item \href{#method-Dataset-api_report}{\code{Dataset$api_report()}}
\item \href{#method-Dataset-api_gwas_release}{\code{Dataset$api_gwas_release()}}
\item \href{#method-Dataset-clone}{\code{Dataset$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-new"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-new}{}}}
\subsection{Method \code{new()}}{
Initialise
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$new(filename = NULL, wd = tempdir(), igd_id = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{filename}}{Path to raw GWAS summary data file}

\item{\code{wd}}{Path to directory to use as the working directory. Will be deleted upon completion - best to keep as the default randomly generated temporary directory}

\item{\code{igd_id}}{Option to provide a specified ID for upload. If none provided then will use the next ieu-a batch ID}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
new ObtainEbiDataset object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-is_new_id"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-is_new_id}{}}}
\subsection{Method \code{is_new_id()}}{
Check if the specified ID is unique within the database. It checks published GWASs and those currently being processed
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$is_new_id(id = self$igd_id)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID to check}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-delete_wd"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-delete_wd}{}}}
\subsection{Method \code{delete_wd()}}{
Delete working directory
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$delete_wd()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-set_wd"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-set_wd}{}}}
\subsection{Method \code{set_wd()}}{
Set working directory (creates)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$set_wd(wd)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{wd}}{working directory}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-se_from_bp"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-se_from_bp}{}}}
\subsection{Method \code{se_from_bp()}}{
Estimate standard error from beta and p-value
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$se_from_bp(beta, pval, minp = 1e-300)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{beta}}{Effect size}

\item{\code{pval}}{p-value}

\item{\code{minp}}{Minimum p-value cutoff default = 1e-300}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-determine_columns"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-determine_columns}{}}}
\subsection{Method \code{determine_columns()}}{
Specify which columns in the dataset correspond to which fields.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$determine_columns(params, nrows = 100, gwas_file = self$filename, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{params}}{List of column identifiers. Identifiers can be numeric position or column header name. Required columns are: c("chr_col", "pos_col", "ea_col", "oa_col", "beta_col", "se_col", "pval_col","rsid_col"). Optional columns are: c("snp_col", "eaf_col", "oaf_col", "ncase_col", "imp_z_col", "imp_info_col", "ncontrol_col").}

\item{\code{nrows}}{How many rows to read to check that parameters have been specified correctly}

\item{\code{gwas_file}}{Filename to read}

\item{\code{...}}{Further arguments to pass to data.table::fread in order to correctly read the dataset}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-format_dataset"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-format_dataset}{}}}
\subsection{Method \code{format_dataset()}}{
Process dataset ready for uploading. Determins build and lifts over to hg19/b37 if necessary.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$format_dataset(
  gwas_file = self$filename,
  gwas_out = file.path(self$wd, "format.txt.gz"),
  params = self$params,
  metadata_test = self$metadata_test,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{gwas_file}}{GWAS filename}

\item{\code{gwas_out}}{Filename to save processed dataset to}

\item{\code{params}}{Column specifications (see determine_columns for more info)}

\item{\code{metadata_test}}{List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats}

\item{\code{...}}{Further arguments to pass to data.table::fread in order to correctly read the dataset}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-view_metadata_options"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-view_metadata_options}{}}}
\subsection{Method \code{view_metadata_options()}}{
View the specifications for available meta data fields, as taken from http://gwas-api.mrcieu.ac.uk/docs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$view_metadata_options()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-get_gwasdata_fields"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-get_gwasdata_fields}{}}}
\subsection{Method \code{get_gwasdata_fields()}}{
Get a list of GWAS data fields and whether or not they are required
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$get_gwasdata_fields()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-get_metadata_fields"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-get_metadata_fields}{}}}
\subsection{Method \code{get_metadata_fields()}}{
Get a list of metadata fields and whether or not they are required
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$get_metadata_fields()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-collect_metadata"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-collect_metadata}{}}}
\subsection{Method \code{collect_metadata()}}{
Input metadata
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$collect_metadata(metadata, igd_id = self$igd_id)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{metadata}}{List of meta-data fields and their values, see view_metadata_options for which fields need to be inputted.}

\item{\code{igd_id}}{ID to be used for uploading to the database}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-check_meta_data"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-check_meta_data}{}}}
\subsection{Method \code{check_meta_data()}}{
Check that the reported effect allele and effect allele frequency columns are correct.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$check_meta_data(
  gwas_file = self$filename,
  params = self$params,
  metadata = self$metadata
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{gwas_file}}{Filename to read}

\item{\code{params}}{column names from x$determine_columns(). Required columns are: c("snp_col", "ea_col", "oa_col", "eaf_col" )}

\item{\code{metadata}}{metadata from x$collect_metadata()}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-write_metadata"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-write_metadata}{}}}
\subsection{Method \code{write_metadata()}}{
Write meta data to json file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$write_metadata(
  metadata = self$metadata,
  datainfo = self$datainfo,
  outdir = self$wd
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{metadata}}{List of meta data fields and their values}

\item{\code{datainfo}}{List of data column parameters}

\item{\code{outdir}}{Output directory to write json files}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_metadata_upload"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_metadata_upload}{}}}
\subsection{Method \code{api_metadata_upload()}}{
Upload meta data to API
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_metadata_upload(
  metadata = self$metadata,
  metadata_test = self$metadata_test,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{metadata}}{List of meta data fields and their values}

\item{\code{metadata_test}}{List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_metadata_edit"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_metadata_edit}{}}}
\subsection{Method \code{api_metadata_edit()}}{
Upload meta data to API
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_metadata_edit(
  metadata = self$metadata,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{metadata}}{List of meta data fields and their values}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_metadata_check"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_metadata_check}{}}}
\subsection{Method \code{api_metadata_check()}}{
View meta-data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_metadata_check(
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID to check}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_metadata_delete"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_metadata_delete}{}}}
\subsection{Method \code{api_metadata_delete()}}{
Delete a dataset. This deletes the metadata AND any uploaded GWAS data (and related processing files)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_metadata_delete(
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID to delete}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_gwasdata_upload"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_gwasdata_upload}{}}}
\subsection{Method \code{api_gwasdata_upload()}}{
Upload gwas dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_gwasdata_upload(
  datainfo = self$datainfo,
  gwasfile = self$gwas_out,
  metadata_test = self$metadata_test,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datainfo}}{List of data column parameters}

\item{\code{gwasfile}}{Path to processed gwasfile}

\item{\code{metadata_test}}{List of outputs from tests of the effect allele, effect allele frequency columns and summary data using CheckSumStats}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_gwasdata_check"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_gwasdata_check}{}}}
\subsection{Method \code{api_gwasdata_check()}}{
Check status of API processing pipeline
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_gwasdata_check(
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID to check}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_gwasdata_delete"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_gwasdata_delete}{}}}
\subsection{Method \code{api_gwasdata_delete()}}{
Delete a dataset. This deletes the metadata AND any uploaded GWAS data (and related processing files)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_gwasdata_delete(
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID to delete}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_qc_status"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_qc_status}{}}}
\subsection{Method \code{api_qc_status()}}{
Check the status of the GWAS QC processing pipeline
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_qc_status(
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID to delete}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_report"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_report}{}}}
\subsection{Method \code{api_report()}}{
View the html report for a processed dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_report(
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{ID of report to view}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-api_gwas_release"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-api_gwas_release}{}}}
\subsection{Method \code{api_gwas_release()}}{
Release a dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$api_gwas_release(
  comments = NULL,
  passed_qc = "True",
  id = self$igd_id,
  access_token = ieugwasr::check_access_token()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{comments}}{Optional comments to provide when uploading}

\item{\code{passed_qc}}{True or False}

\item{\code{id}}{ID to release}

\item{\code{access_token}}{Google OAuth2.0 token. See ieugwasr documentation for more info}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Dataset-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Dataset-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Dataset$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
