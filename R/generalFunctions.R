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
    if (length(names(l)) != length(l)) {
        stop("The list you submitted might not be named.")
    }
    if (!is.null(names(l[[1]]))) {
        ni <- gsub(".*\\.(.*)$", "\\1", names(unlist(l)))
        n <- rep(names(l), lengths(l))
        tibble(
            names = n,
            inner.names = ni,
            variables = unname(unlist(l))
        )
    } else {
        n <- rep(names(l), lengths(l))
        tibble(
            names = n,
            variables = unname(unlist(l))
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
#' @importFrom AnnotationDbi select keys
#' @import GO.db

getGO.db <- function() {
  AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
}

#' col64
#'
#' Diverging color palette.
#'
#' @name col64
#' @rdname col64
#' @author Jason T. Serviss
#' @keywords col64
#' @examples
#' col64()
#'
#' @export


col64 <- function() {
    c(
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6",
    "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43",
    "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601",
    "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0",
    "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B",
    "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500",
    "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C",
    "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
    "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
    "#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700",
    "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837",
    "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625",
    "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC",
    "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"
    )
}

#' Get ggplot legend.
#'
#'
#' @name g_legend
#' @rdname g_legend
#' @author Jason T. Serviss
#' @keywords g_legend
#' @param a.gplot A ggplot.
#' @examples
#' #no example yet
#'
#' @export
#' @importFrom ggplot2 ggplot_gtable ggplot_build

g_legend <- function(a.ggplot){
    tmp <- ggplot_gtable(ggplot_build(a.ggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
