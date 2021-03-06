#' Gene ontology term enrichment results.
#'
#' @title Results of the gene ontology term enrichment analysis with topGO.
#' @description Produced with \code{\link{runTopGO}}.
#' @docType data
#' @name topGOresult
#' @format data frame
#' \describe{
#'   \item{GO.ID}{character. Gene ontology term ID.}
#'   \item{Term}{character. Gene ontology term description.}
#'   \item{Annotated}{int. Derived from \link[topGO]{GenTable}.}
#'   \item{Significant}{int. Derived from \link[topGO]{GenTable}.}
#'   \item{Expected}{numeric. Derived from \link[topGO]{GenTable}.}
#'   \item{classicFisher}{numeric. Derived from \link[topGO]{GenTable}.}
#'   \item{gene_name}{character. A comma seperated vector of gene names related
#'        to the GO term.}
#'   \item{ID}{character. A comma seperated vector of Ensembl IDs related
#'        to the GO term.}
#'   \item{biotype}{character. A comma seperated vector of gene biotypes related
#'        to the GO terms genes.}
#' }
#' @usage topGOresult
#' @return A data frame containing the results from the
#'    \href{http://bioconductor.org/packages/topGO/}{topGO} gene ontology term
#'    enrichment analysis. The first 6 columns are derived from the
#'    \link[topGO]{GenTable} function. Genes related to each GO term were
#'    annotated with the \link[topGO]{genesInTerm} function and further gene
#'    information was merged from the \code{\link{geneAnnotation}} dataset.
#' @examples
#' topGOresult
#'
NULL
