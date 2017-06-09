#!/usr/bin/env Rscript
#' Run DESeq2 and limma for HCT116 cell line
#'
#'
#' This command runs DESeq2 using output by HTSeq for the HCT116 experiment.
#'    It compares parental and adapted cell lines and estimates the differential
#'    expression between the two using the DESeq2 software.
#'    Furthermore, the results are annotated with normalized expression counts,
#'    and gene annotation.
#'
#' @name runHCT
#' @rdname runHCT
#' @author Jason T. Serviss
#' @keywords runHCT
#' @examples
#'
#' \dontrun{runHCT()}
NULL

runHCT <- function() {
    
    packagePath <- getwd()
  
    # directories, paths, source functions
    dataINdir <- file.path(packagePath,'inst', 'data-raw')
    dataOUTdir <- file.path(packagePath, 'data')
    source(file.path(dataINdir, 'deFunctions.R'))
  
    ##check for annotation and load
    gtfSimplePATH <- file.path(packagePath, 'data', 'annotation.rda')
  
    if(! file.exists(gtfSimplePATH)){
        print("making annotation")
        c <- "Rscript -e \"source('./inst/R.scripts/getAnnotation.R'); getAnnotation()\""
        system(c)
    }
  
    load(gtfSimplePATH)
  
    # prepare the counts and samples info
    #sample info
    library_info <- data.frame(
        sample <- c(
            "HCT116_pH7A",
            "HCT116_pH7B",
            "HCT116_pH7C",
            "HCT116_pH6A",
            "HCT116_pH6B",
            "HCT116_pH6C"
        ),
        library <- c(
            "parentalA",
            "parentalB",
            "parentalC",
            "adaptedA",
            "adaptedB",
            "adaptedC"
        ),
        condition <- factor(
            c(
                rep("parental", 3),
                rep("adapted", 3)
            ),
            levels=c("parental", "adapted")
        )
    )
    targets <- library_info[,(2:3)]
    SLLids <- as.vector(library_info[,1])
  
    #counts
    input_file <- file.path(dataINdir, "count_table.txt")
    counts <- read.table(input_file, head=T, sep="\t")
    counts <- counts[, c(
        "HCT116_pH7A", "HCT116_pH7B", "HCT116_pH7C", "HCT116_pH6A",
        "HCT116_pH6B", "HCT116_pH6C", "Me3pH7A", "Me3pH7B", "Me3pH7C",
        "Me3pH6A", "Me3pH6B", "Me3pH6C")]
    counts <- counts[, 1:6]
  
  
    counts <- counts[,SLLids]
    libraries <- targets$library
    colnames(counts) <- libraries
    condition <- targets$condition
    samples <- levels(targets$condition)
    exp <- as.data.frame(targets[,2])
    rownames(exp) <- targets[,1]
    colnames(exp) <- "condition"
  
    ##filter expressed counts
    counts_filtered <- .filter_expressed(counts,1,3)
  
    ##############################################
    #
    #DE using DESeq2
    #
    ##############################################
  
    print("running DESeq2")
    exp.cnts <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_filtered,
        colData = exp,
        design = ~ condition
    )
    exp.deg <- DESeq2::DESeq(exp.cnts) ## fitType="parametric"
    results <- DESeq2::results(exp.deg)
    results$ID <- rownames(results)
    DESeq2 <- as.data.frame(results[,c(7,1:6)])
  
    ##save rld
    rld <- DESeq2::rlogTransformation(exp.deg)
    DESeq2rld.HCT116 <- as.data.frame(SummarizedExperiment::assay(rld))
    colnames(DESeq2rld.HCT116) <- paste(
        "HCT116",
        colnames(DESeq2rld.HCT116),
        sep="."
    )
    fileOUT4 <- paste(dataOUTdir, "deResultsRld.rda", sep="/")
    deResultsRld <- DESeq2rld.HCT116
    save(deResultsRld, file=fileOUT4)
  
    ##merge with counts and gtf
    norm.counts <- as.data.frame(
        DESeq2::counts(exp.deg,
        normalized=TRUE)
    )
    norm.counts$ID <- rownames(norm.counts)
    DESeq2 <- merge(DESeq2, norm.counts, by = "ID", all.x=TRUE)
    DESeq2 <- merge(DESeq2, annotation, by = "ID", all.x=TRUE)
    DESeq2 <- DESeq2[, c(
        "ID", "chr", "start", "end", "geneName", "biotype", "parentalA",
        "parentalB", "parentalC", "adaptedA", "adaptedB", "adaptedC",
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    DESeq2.HCT116 <- DESeq2[order(DESeq2$padj), ]
  
    ##save
    hsa.res.pth <- file.path(dataOUTdir, 'deResults.rda')
    deResults <- DESeq2.HCT116
    save(deResults, file=hsa.res.pth)
  
    print("done")
}







