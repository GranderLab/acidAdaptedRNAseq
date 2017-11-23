
#' TSNE
#'
#' Runs t-SNE and plots the results.
#'
#' This function runs t-SNE for all samples and the specified genes using
#' 1-Pearson's correlation as input to the t-SNE algorithm and plots the results.
#'
#' @name TSNE
#' @rdname TSNE
#' @author Jason T. Serviss
#' @param ntop The number of genes to include in t-SNE ranked by their variance.
#' @examples
#' TSNE()
#'
NULL

#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename mutate
#' @importFrom stringr str_sub
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_ptol
#' @importFrom stats var
#' @importFrom magrittr "%>%"

TSNE <- function(nTop = 500) {
  x <- deResultsRld
  
  #select top genes
  rv <- apply(x, 1, var)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  #calculate 1-Pearson's correlation
  d <- as.dist(1-cor(x[select, ], method = "p"))
  
  #run tsne, fix labels, and plot
  as.data.frame(Rtsne(d, is_distance = TRUE, perplexity = 1)$Y) %>%
  as_tibble() %>%
  rename(dim1 = V1, dim2 = V2) %>%
  mutate(samples = gsub("\\.", " ", colnames(x))) %>%
  mutate(samples = gsub("(.*)(.)$", "\\1 \\2", samples)) %>%
  mutate(samples = case_when(
    str_sub(samples, start = -1) == "A" ~ gsub("(.*).$", "\\1rep 1", samples),
    str_sub(samples, start = -1) == "B" ~ gsub("(.*).$", "\\1rep 2", samples),
    str_sub(samples, start = -1) == "C" ~ gsub("(.*).$", "\\1rep 3", samples),
    TRUE ~ "error"
  )) %>%
  mutate(group = gsub("(.*).{6}$", "\\1", samples)) %>%
  ggplot(aes_string('dim1', 'dim2', 'colour = group')) +
    geom_point(size = 3, alpha = 0.7) +
    theme_few() +
    scale_colour_ptol() +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title = "Condition"))
}

#' MAplot
#'
#' Plot an MA plot (log 2 average expression vs. log2 fold change) using the
#' results from the DESeq2 differential expression analysis.
#'
#' @name MAplot
#' @rdname MAplot
#' @author Jason T. Serviss
#' @keywords MAplot
#' @examples
#'
#' MAplot()
#'
#' @export
#' @importFrom tibble as_tibble
#' @importFrom dplyr select
#' @import ggplot2
#' @importFrom ggthemes theme_few
#' @importFrom magrittr "%>%"
NULL

MAplot <- function() {
    data <- deResults %>%
        as_tibble() %>%
        select(~baseMean, ~log2FoldChange)
    
    p <- ggplot(data, aes_string('log2(baseMean)', 'log2FoldChange'))+
        geom_point(alpha = 0.1)+
        theme_few()+
        labs(
            x = "log2(base mean)",
            y = "log2(fold change)",
            title = "MA plot",
            subtitle = "Red lines correspond to a fold change of 2."
        )+
        ylim(-6, 6)+
        geom_hline(yintercept = 1, colour = "red", lty = 2)+
        geom_hline(yintercept = -1, colour = "red", lty = 2)+
        theme(
            plot.subtitle = element_text(size = 10)
        )
    
    p
    return(p)

}

#' volcanoPlot
#'
#' Plot an volcano plot using the results from the DESeq2 differential
#' expression analysis using only significant (alpha < 0.05) genes.
#'
#' @name volcanoPlot
#' @rdname volcanoPlot
#' @author Jason T. Serviss
#' @param n Numeric. The number of gene names to plot based on their log2 fold
#'    change.
#' @keywords volcanoPlot
#' @examples
#'
#' volcanoPlot()
#'
#' @export
#' @importFrom tibble as_tibble
#' @importFrom dplyr select_ filter_ bind_rows arrange_ desc
#' @import ggplot2
#' @importFrom ggthemes theme_few
#' @importFrom utils head
#' @importFrom magrittr %>%
NULL

volcanoPlot <- function(n = 5) {
    data <- deResults %>%
        as_tibble() %>%
        filter_(~padj < 0.05) %>%
        select_(~log2FoldChange, ~padj)
    
    genesTop <- deResults %>%
        as_tibble() %>%
        filter_(~padj < 0.05) %>%
        arrange_(~desc(log2FoldChange)) %>%
        head(n = n) %>%
        select_(~log2FoldChange, ~padj, ~geneName)
    
    genesBottom <- deResults %>%
        as_tibble() %>%
        filter_(~padj < 0.05) %>%
        arrange_(~log2FoldChange) %>%
        head(n = n) %>%
        select_(~log2FoldChange, ~padj, ~geneName)
    
    genes <- bind_rows(genesTop, genesBottom)
    
    p <- ggplot(data, aes_string('log2FoldChange', '-log10(padj)'))+
        geom_point(alpha = 0.1)+
        geom_text(
            data = genes,
            aes_string(
                x = 'log2FoldChange', 
                y = '-log10(padj)', 
                label = 'geneName'
            ),
            nudge_y = 1
        )+
        theme_few()+
        geom_vline(xintercept = -1, colour = "red", lty = 2)+
        geom_vline(xintercept = 1, colour = "red", lty = 2)+
        labs(
            x = "log2(fold change)",
            y = "-log10(p-value)",
            title = "Volcano plot",
            subtitle = "Red lines represent a fold change of 2."
        )+
        theme(
            plot.subtitle = element_text(size = 10)
        )+
        xlim(
            max(range(data$log2FoldChange)) * -1,
            max(range(data$log2FoldChange))
        )
    p
    return(p)
}
