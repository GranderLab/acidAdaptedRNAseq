
getAnnotation <- function() {
  tmp <- tempdir()
  
    sCmd1 <- paste(
        "wget ",
        "ftp://ftp.ensembl.org/",
        "pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz ",
        "-O ", tmp, "/Homo_sapiens.GRCh37.73.gtf.gz",
        sep = ""
    )
    system(sCmd1)
    
    sysCmd2 <- paste(
      "gunzip -f ",
      file.path(tmp, "Homo_sapiens.GRCh37.73.gtf.gz"),
      sep = "")
    system(sysCmd2)
    
    library(GenomicFeatures)
    txdb <- makeTxDbFromGFF(
        file.path(tmp, "Homo_sapiens.GRCh37.73.gtf"),
        format = "gtf"
    )
    annotation <- as.data.frame(genes(txdb))[, -4]
    colnames(annotation) <- c('chr', 'start', 'end', 'strand', 'ID')
    
    otherData <- .getOtherData(annotation)
    annotation <- merge(annotation, otherData, by = "ID")
    
    set <- c(
        'ID',
        'chr',
        'start',
        'end',
        'strand',
        'geneName',
        'biotype',
        'description'
    )
    
    annotation <- annotation[, set]
    save(annotation, file = "./data/annotation.rda", compress = "bzip2")
}

.getOtherData <- function(annotation){
    values <- annotation$ID
    ensembl <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        #host = "sep2013.archive.ensembl.org"
        host = "dec2013.archive.ensembl.org"
    )
    
    IDs <- biomaRt::getBM(
        attributes = c(
            'ensembl_gene_id',
            'description',
            'gene_biotype',
            'external_gene_id'
        ),
        filters = 'ensembl_gene_id',
        values = values,
        mart = ensembl
    )
    
    colnames(IDs) <- c("ID", "description", "biotype", "geneName")
    IDs$description <- gsub("(.*) \\[Source:.*\\]", "\\1", IDs$description)
    return(IDs)
}
