#' runTopGO
#'
#' Runs GO analysis using topGO.
#'
#' @name runTopGO
#' @rdname runTopGO
#' @author Jason T. Serviss
#' @param save Logical. Indicates if the results should be
#'    saved to ./data/topGOresult.rda.
#' @keywords runTopGO
#' @examples
#' \dontrun{runTopGO()}
#'
#' @export
#' @import topGO
#' @importFrom biomaRt useMart getBM
#' @importFrom methods new
NULL

runTopGO <- function(save=FALSE) {
    
    ##subset gene universe which is comprised of all quantified genes
    universe <- unique(deResults$ID)
    
    #only use DE in HCT116 found with DESeq2
    DE <- deResults[deResults$padj < 0.05, "ID"]
    
    ##runGO
    GO = .runGO(DE, universe)
    
    ##get full GO names
    ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    go <- getBM(attributes=c('go_id', 'name_1006'), mart = ensembl)
    GO <- merge(GO, go, by.x="GO.ID", by.y="go_id", all.x=TRUE, all.y=FALSE)
    GO$name_1006[is.na(GO$name_1006)] <- GO$Term[is.na(GO$name_1006)]
    GO$Term <- GO$name_1006
    GO$name_1006 <- NULL
    
    ##check coverage
    .Pcoverage(GO, DE)
    
    if(save) {
        save(GO, file="data/topGOresult.rda", compress = "bzip2")
    }
    
    return(GO)
}

.runGO <- function(DE, universe) {
    
    ##reformat DE results for topGO
    geneList <- factor(as.integer(universe %in% DE))
    names(geneList) <- universe
    
    ##run GO
    GOdata <- new(
        "topGOdata",
        ontology = "BP",
        allGenes = geneList,
        nodeSize = 5,
        annot = annFUN.org,
        mapping = "org.Hs.eg.db",
        ID = "Ensembl"
    )
    
    resultFisher <- runTest(
        GOdata,
        algorithm = "classic",
        statistic = "fisher"
    )
    
    outTmp <- GenTable(
        GOdata,
        classicFisher = resultFisher,
        topNodes = 1000
    )
    
    outTmp2 <- outTmp[outTmp$classicFisher < 1, ]
    
    ##annotate genes associated with reported GO terms
    allGO = genesInTerm(GOdata) ###gets all genes related to go terms
    newGO = .namedListToDF(x = allGO, y = DE)
    
    ##merge GO results to include the genes involved in the GO terms
    merg = merge(outTmp2, newGO, by = "GO.ID", all.x = TRUE)
    merg = merg[order(merg$classicFisher), ]
    return(merg)
}

##for use with the runGO function.
#Turns a named list of unequal length into a data frame.
#This also annotates gene name from given gene ID.

.namedListToDF <- function(x, y) {
    out = data.frame()
    
    for( tt in 1:length(x) ) {
        merg <- data.frame()
        curr = x[tt]
        name = names(curr)
        INvalues = data.frame(ID = unlist(strsplit(curr[[1]], ", ")))
        merg = merge(geneAnnotation, INvalues, by="ID", all.y=TRUE)
        subs <- merg[merg$ID %in% y, ]
        out[tt, "GO.ID"] <- name
        out[tt, "gene_name"] <- ifelse(
            nrow(subs) == 0,
            NA,
            paste(subs$geneName, collapse = ", ")
        )
        out[tt, "ID"] <- ifelse(
            nrow(subs) == 0,
            NA,
            paste(subs$ID, collapse = ", ")
        )
        out[tt, "biotype"] <- ifelse(
            nrow(subs) == 0,
            NA,
            paste(subs$biotype, collapse = ", ")
        )
    }
    return(out)
}

##calculate coverage. for use with the runGO function.
.Pcoverage <- function(x, interesting) {
    
    output <- data.frame()
    
    ##calculate number of significant GO terms
    sig.nr <- nrow(x[x$classicFisher < 0.05, ])
    
    ##calculate most significant GO term
    most.sig <- min(x$classicFisher)
    return <- unique(unlist(strsplit(x$ID, ", ")))
    
    ##calculate coverage and print output
    final <- unique(return)
    found <- final[final %in% interesting]
    percent <- (length(final) / length(interesting))*100
    
    print("% coverage:")
    print(percent)
    print("total DE genes:")
    print(length(interesting))
    print("genes in GO:")
    print(length(final))
    print("nr. significant GO terms:")
    print(sig.nr)
    print("most significant term:")
    print(most.sig)
    cat("\n")
    
    return(output)
}

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
NULL

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
#' @param spins The maximum number of spins allowed. Passed
#'    to \link[igraph]{cluster_spinglass}.
#' @keywords collapseGraph
#' @return A list with the first element being the collapsed graph and the
#'    second element being a summary of the community analysis.
#' @examples
#' terms <- c("GO:0006915", "GO:0097152", "GO:0006664", "GO:1903509")
#' g <- downloadGOgraph(terms)
#' collapseGraph(g)
#'
#' @export
#' @importFrom igraph cluster_spinglass authority_score contract.vertices V
#'    membership simplify
#' @importFrom tibble tibble rownames_to_column add_column
#' @importFrom stats setNames
#' @importFrom dplyr left_join pull select_ mutate n
#' @importFrom magrittr %>%
NULL

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
        mapping = m,
        vertex.attr.comb=function(x) maxAuth(x, as)
    ) %>%
    igraph::simplify(edge.attr.comb=list(weight="sum"))
}

.assembleCommunitySummary <- function(spinglass, as, cg) {
    
    #get GO annotaion
    vals <- getGO.db()
    
    #process authority score
    as <- as %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      as_tibble() %>%
      setNames(c("vertex", "score"))
    
    #process spinglass
    l <- lapply(1:length(spinglass), function(x) spinglass[[x]])
    names(l) <- as.character(1:length(l))
    
    #process community IDs
    comID <- V(cg)$name %>%
        as_tibble() %>%
        setNames("communityID") %>%
        mutate(communityName = as.character(1:n()))
        
    #assemble summary
    summaryNames <- c(
        "communityName",
        "GOID",
        "authority.score",
        "communityID",
        "Term",
        "communityTerm",
        "communitySize"
    )
    
    . <- NULL
    summary <- namedListToTibble(l) %>%
      left_join(as, by=c("variables" = "vertex")) %>%
      left_join(comID, by=c("names" = "communityName")) %>%
      left_join(vals, by=c("variables" = "GOID")) %>%
      select_(~-ONTOLOGY) %>%
      left_join(vals, by=c("communityID" = "GOID")) %>%
      select_(~-ONTOLOGY) %>%
      add_column(size = spinglass$csize[as.numeric(pull(., names))]) %>%
      setNames(summaryNames) %>%
      select_(
        ~communityName,
        ~communityID,
        ~communitySize,
        ~communityTerm,
        ~GOID,
        ~authority.score,
        ~Term
      )
    
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
#' @param assignment A data frame of edges from term to community term with
#'    colnames \emph{from} and \emph{to}.
#' @keywords addBackVertexAndEdges
#' @examples
#' terms <- c("GO:0006915", "GO:0097152", "GO:0006664", "GO:1903509")
#' g <- downloadGOgraph(terms)
#' collapseGraph(g)
#'
#' @export
#' @importFrom igraph get.vertex.attribute add_vertices add_edges
#' @importFrom GOSim getGOGraph
#' @importFrom magrittr %>%
NULL

addBackVertexAndEdges <- function(g, cg, assignment) {
    bools1 <- !get.vertex.attribute(g)$name %in% get.vertex.attribute(cg)$name
    vertexToAdd <- get.vertex.attribute(g)$name[bools1]
    
    cg <- cg %>%
        add_vertices(length(vertexToAdd), attr = list(name=vertexToAdd)) %>%
        add_edges(c(rbind(
            match(assignment$from, get.vertex.attribute(.)$name),
            match(assignment$to, get.vertex.attribute(.)$name)
            )
        ))

}






