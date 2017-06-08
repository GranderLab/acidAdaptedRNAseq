#' downloadGOgraph
#'
#' Downloads the GO graph for a set of GO IDs and converts it to igraph format.
#'
#' @name downloadGOgraph
#' @rdname downloadGOgraph
#' @author Jason T. Serviss
#' @param terms Character. A character vector of GO IDs.
#' @keywords downloadGOgraph
#' @return A GO graph in igraph format with vertex attributes \emph{name} set to
#'    the GO IDs and the edge attribute \emph{weight} set to 1 for all edges.
#' @examples
#' downloadGOgraph("GO:0006915")
#'
#' @export
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom GOSim getGOGraph

downloadGOgraph <- function(terms) {
    igraph.from.graphNEL(getGOGraph(terms))
}

#' collapseGraph
#'
#' Takes a GO graph in igraph format, runs the spin-glass community algorithm,
#' and collapses the graph by community using the name of the node with the
#' highest authority score to name each vertex in the community.
#'
#' @name collapseGraph
#' @rdname collapseGraph
#' @author Jason T. Serviss
#' @param g A GO graph in igraph format.
#' @keywords collapseGraph
#' @return A list with the first element being the collapsed graph and the
#'    second element being a summary of the community analysis.
#' @examples
#' terms <- c("GO:0006915", "GO:0097152", "GO:0006664", "GO:1903509")
#' g <- downloadGOgraph(terms)
#' collapseGraph(g)
#'
#' @export
#' @importFrom igraph cluster_spinglass authority_score contract.vertices
#'    membership simplify
#' @importFrom tibble tibble

collapseGraph <- function(g, spins=200) {
    
    #run community analysis
    set.seed(1998)
    spinglass <- cluster_spinglass(g, spins=200)
    
    #calculate the authority score for each vertex in the graph
    as <- authority_score(g)$vector
    
    #collapse and simplify the graph and name each new vertex according to the
    #vertex with the highest authority score per community.
    cg <- .contractAndSimplify(g, spinglass, as)
    
    #assemble summary
    summary <- .assembleCommunitySummary(spinglass, as)
    
    return(list(cg, summary))
}

#determines the vertex with the highest authority score per community
maxAuth <- function(x, as) {
    curr <- as[names(as) %in% x]
    names(curr)[curr == max(curr)][1]
}

#collapse and simplify the graph
.contractAndSimplify <- function(g, spinglass, as) {
    m <- membership(spinglass)
    
    cg <- g %>%
      contract.vertices(
        graph = .,
        mapping = m,
        vertex.attr.comb=function(x) maxAuth(x, as)
    ) %>%
    igraph::simplify(., edge.attr.comb=list(weight="sum"))
}

.assembleCommunitySummary <- function(spinglass, as) {
    
    #process authority score
    as <- as %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      as_tibble() %>%
      setNames(., c("vertex", "score"))
    
    #process spinglass
    l <- lapply(1:length(spinglass), function(x) spinglass[[x]])
    names(l) <- as.character(1:length(l))
    
    #assemble summary
    summaryNames <- c(
        "communityName",
        "GOID",
        "authority.score",
        "communitySize"
    )
    
    summary <- namedListToTibble(l) %>%
      left_join(as, by=c("variables" = "vertex")) %>%
      add_column(size = spinglass$csize[as.numeric(pull(., names))]) %>%
      setNames(., summaryNames)
    
    return(summary)
}

#' addBackVertexAndEdges
#'
#' Takes a GO graph in igraph format, runs the spin-glass community algorithm,
#' and collapses the graph by community using the name of the node with the
#' highest authority score to name each vertex in the community.
#'
#' @name addBackVertexAndEdges
#' @rdname addBackVertexAndEdges
#' @author Jason T. Serviss
#' @param g A GO graph in igraph format.
#' @param cg The collapsed graph from the function \code{\link{collapseGraph}}.
#' @param assignment A data frame of edges from term to community term.
#' @keywords addBackVertexAndEdges
#' @examples
#' terms <- c("GO:0006915", "GO:0097152", "GO:0006664", "GO:1903509")
#' g <- downloadGOgraph(terms)
#' collapseGraph(g)
#'
#' @export
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom GOSim getGOGraph

addBackVertexAndEdges <- function(g, cg, assignment) {
    bools1 <- !get.vertex.attribute(g)$name %in% get.vertex.attribute(cg)$name
    vertexToAdd <- get.vertex.attribute(g)$name[bools1]
    
    cg <- cg %>%
    add_vertices(., length(vertexToAdd), attr = list(name=vertexToAdd)) %>%
    add_edges(., c(rbind(match(edgesToAdd$GOID, get.vertex.attribute(.)$name), match(edgesToAdd$groupGOID, get.vertex.attribute(.)$name))))

}






