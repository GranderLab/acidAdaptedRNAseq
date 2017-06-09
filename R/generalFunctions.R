#' namedListToTibble
#'
#' Converts a named list to a long data frame.
#'
#' @name namedListToTibble
#' @rdname namedListToTibble
#' @author Jason T. Serviss
#' @param l List. The list to be converted.
#' @keywords namedListToTibble
#' @examples
#'
#' l <- list(a=LETTERS[1:10], b=letters[1:5])
#' namedListToTibble(l)
#'
#' @export
#' @importFrom tibble tibble

namedListToTibble <- function(l) {
    if(length(names(l)) != length(l)) {
        stop("The list you submitted might not be named.")
    }
    if(!is.null(names(l[[1]]))) {
        ni <- gsub(".*\\.(.*)$", "\\1", names(unlist(l)))
        n <- rep(names(l), lengths(l))
        tibble::tibble(
            names=n,
            inner.names=ni,
            variables=unname(unlist(l))
        )
    } else {
        n <- rep(names(l), lengths(l))
        tibble::tibble(
            names=n,
            variables=unname(unlist(l))
        )
    }
}

#' getGO.db
#'
#' Get annotation from \link[GO.db]{GO.db}.
#'
#' @name getGO.db
#' @rdname getGO.db
#' @author Jason T. Serviss
#' @keywords getGO.db
#' @examples
#' getGO.db()
#'
#' @export
#' @importFrom AnnotationDbi select
#' @import GO.db

getGO.db <- function() {
    select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
}