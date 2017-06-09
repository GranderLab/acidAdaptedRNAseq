#' Gene Annotation.
#'
#' @title Gene nnotation from ensembl.
#' @description Assembled with data-raw/getAnnotation.R getAnnotation().
#' @docType data
#' @name geneAnnotation
#' @format data frame
#' \describe{
#'  \item{ID}{character. Ensembl ID.}
#'  \item{chr}{character. Chromosome.}
#'  \item{start}{int. Transcription start site.}
#'  \item{end}{int. Transcription end site.}
#'  \item{strand}{factor. Strand.}
#'  \item{geneName}{character. Gene name.}
#'  \item{biotype}{character. Gene biotype.}
#'  \item{description}{character. Gene description.}
#' }
#' @usage geneAnnotation
#' @return A data frame with Ensembl gene annotation.
#' @examples
#' geneAnnotation
#'

NULL