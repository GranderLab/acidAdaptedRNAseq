#' Regularized log transformed expression values.
#'
#' @title Results from the DESeq2 differential expression analysis.
#' @description Produced with data-raw/runHCT.R runHCT().
#' @docType data
#' @name deResultsRld
#' @format data.frame
#' \describe{
#'   \item{rownames}{Ensembl IDs}
#'   \item{colnames}{Sample identification}
#' }
#' @usage deResultsRld
#' @return A data frame containing the regularized log transformed counts after
#'    differential expression analysis with
#'    \href{http://bioconductor.org/packages/DESeq2/}{DESeq2}. The data was
#'    generated using the \link[DESeq2]{rlogTransformation} function.
#' @examples
#' deResultsRld
#'
NULL