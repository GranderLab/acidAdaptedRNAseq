#'
#' Functions for qc of RNAseq data and detection of differentially expressed
#' genes.
#'
#' Internal function for differential expression analysis
#'
#' filter out features with less than x cpm count in more than n libraries
#' default x = 1 cpm and n is the number in smallest replicate.
#' recommended by Gordon Smyth (edgeR)
#' counts is a df with feature IDs as rownames
#'
#' @name .filter_expressed
#' @rdname .filter_expressed
#' @author Agata Smialowska
#' @keywords internal
#' @examples
#'
#' ##no example yet
#'


.filter_expressed <- function(counts, x, n) {
    cpm <- as.matrix(counts)
    cpm <- t(t(cpm) / colSums(cpm) * 1000000)
    is.expressed_1 <- rowSums(cpm > x) >= n
    counts_filt <- counts[is.expressed_1, ]
    
    return(counts_filt)
}

#'
#' Functions for qc of RNAseq data and detection of differentially expressed genes
#'
#' Internal function for differential expression analysis
#'
#' contrasts for limma
#' process metadata to contrasts and generate design matrix (for limma)
#' metadata is in df with library & condition columns
#' samples=levels(targets$condition)
#'
#' @name .design.contrasts.limma
#' @rdname .design.contrasts.limma
#' @author Agata Smialowska
#' @keywords internal
#' @examples
#'
#' ##no example yet
#'


.design.contrasts.limma <- function(targets){
    
    samples.cont <- levels(targets$condition)
    f <- factor(targets$condition, levels = samples.cont)
    design <- model.matrix(~0 + f)
    colnames(design) <- samples.cont
    
    return(design)
    
}

#'
#' Functions for qc of RNAseq data and detection of differentially expressed
#' genes
#'
#' Internal function for differential expression analysis
#'
#' contrasts for limma
#' process metadata to contrasts and generate design matrix (for limma)
#' metadata is in df with library & condition columns
#' samples=levels(targets$condition)
#'
#' @name .contrasts.limma
#' @rdname .contrasts.limma
#' @author Agata Smialowska
#' @keywords internal
#' @examples
#'
#' ##no example yet
#'


.contrasts.limma <- function(targets, design){
    
    samples.cont <- levels(targets$condition)
    contrast_group <- list()
    
    for (i in samples.cont) {
        con.1 <- i
        for (i in samples.cont) {
            con.2 <- i
            if (con.2 != con.1) {
                contrast <- paste(con.1, "-", con.2,sep = "")
                contrast_group <- paste(contrast_group, contrast, sep = " ")
            }
        }
    }
    
    contrasts <- strsplit(
        contrast_group,
        split = " ",
        fixed = TRUE,
        perl = FALSE,
        useBytes = FALSE
    )
    contrasts <- unlist(contrasts)
    conts <- grep("-", contrasts, value = TRUE)
    contrast.matrix <- limma::makeContrasts(contrasts = conts, levels = design)
    
    return(contrast.matrix)
    
}