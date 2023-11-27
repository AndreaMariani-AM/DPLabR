#' Volcano Plot
#'
#' @description Custom code that creates a \code{\link[ggplot2]{ggplot}} object plotting differentially
#'            expressed genes.
#' @param df A dataframe containing differentally expressed genes as an output from \code{\link[DESeq2]{DESeq}}.
#'            For the DP lab, this is the output of the contrast section of the pipeline.
#' @param xlim User defined xlim
#' @param ylim User defined ylim
#' @param main User defined main title
#' @param labelSize labelSize of the number or percentage of Upregulated or Downregulated genes
#' @param pval Statistical threshold to draw pvalue significance line. Default is 0.1
#' @param log2FC Statistical threshold to draw log3FoldChange line. Default is 0.1
#' @param raw_numbers Logical. Whether to draw the number of UP/DOWN genes as an integer or a fraction
#'
#' @importFrom rlang .data
#' @import ggplot2
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @export
#'
VolcanoPlot <- function(df, xlim=NULL,
                        ylim=NULL,
                        main = NULL,
                        labelSize = 7,
                        pval = 0.1, log2FC = 1,
                        raw_numbers = TRUE) {

  p <-  ggplot2::ggplot(data = df, aes(x=.data$log2FoldChange, y=-log10(.data$padj), colour = .data$DEG) ) +
    ggplot2::geom_point(alpha=0.7, size=2)

  if (raw_numbers == TRUE) {

    # Get the full number
   p + ggplot2::annotate("text", label = sum(df$DEG == "Upregulated"),
                         color = "red", y = 0, x = xlim[2],
                         vjust="inward",hjust="inward", size = labelSize) +
      ggplot2::annotate("text", label = sum(df$DEG == "Downregulated"),
                        color = "darkgreen", y = 0, x = xlim[1],
                        vjust="inward",hjust="inward", size = labelSize)

  } else {

    # Get the percentage
    p + ggplot2::annotate("text", label = paste(round(sum(df$DEG == "Upregulated")/length(df$DEG)*100,0), "%"),
                          color = "red", y = 0, x = xlim[2],
                          vjust="inward",hjust="inward", size = labelSize) +
      ggplot2::annotate("text", label = paste(round(sum(df$DEG == "Downregulated")/length(df$DEG)*100,0), "%"),
                        color = "darkgreen", y = 0, x = xlim[1],
                        vjust="inward",hjust="inward", size = labelSize)
  }

  p +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.position = "top") +


  ggplot2::ggtitle(main) +
  ggplot2::theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

  ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +

  ggplot2::geom_hline(yintercept = -log10(.data$pval), linetype = 2) +
  ggplot2::geom_vline(xintercept = c(-.data$log2FC, .data$log2FC), linetype = 2) +

  ggplot2::xlab("log2 fold change") + ggplot2::ylab("-log10 p-value") +
  ggplot2::scale_colour_manual(values = c("Downregulated" = "darkgreen", "NS" = "gray", "Upregulated" = "red"),
                        labels = c("Downregulated" = "Dowregulated", "NS" = "NS", "Upregulated" = "Upregulated"),
                        drop = FALSE)

  p
}


#' Volcano Plot with Genes Info
#' @description Creates a \code{\link[ggplot2]{ggplot}} object plotting differentially
#'            expressed genes and highlighting a user defined list of genes
#'
#' @param df A  dataframe containing differentally expressed genes as an output from \code{\link[DESeq2]{DESeq}}.
#'            For the DP lab, this is the output of the contrast section of the pipeline.
#' @param xlim User defined xlim
#' @param ylim User defined ylim
#' @param main User defined main title
#' @param labelSize labelSize of the number or percentage of Upregulated or Downregulated genes
#' @param pval Statistical threshold to draw pvalue significance line. Default is 0.1
#' @param log2FC Statistical threshold to draw log3FoldChange line. Default is 0.1
#' @param Genes A user defined list of genes of interests to be highlighted
#'
#' @importFrom rlang .data
#' @import ggplot2
#' @import ggrepel
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @export
#'

VolcanoPlot_Genes <- function(df,
                              xlim=NULL,
                              ylim=NULL,
                              main = NULL,
                              labelSize = 7,
                              pval = 0.1,
                              log2FC = 1,
                              Genes=NULL) {

  p <-  ggplot2::ggplot(data = df, aes(x=.data$log2FoldChange, y=-log10(.data$padj), colour=.data$DEG) ) +
    ggplot2::geom_point(alpha=0.7, size=2) +

    ggplot2::annotate("text", label = sum(df$DEG == "Upregulated"),
                      color = "red", y = 0, x = xlim[2],
                      vjust="inward",hjust="inward", size = labelSize) +
    ggplot2::annotate("text", label = sum(df$DEG == "Downregulated"),
                      color = "darkgreen", y = 0, x = xlim[1],
                      vjust="inward",hjust="inward", size = labelSize) +

    ggplot2::theme_classic() +
    ggplot2::theme(legend.title = element_blank()) +
    ggplot2::theme(legend.position = "top") +

    ggplot2::ggtitle(main) +
    ggplot2::theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

    ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +

    ggplot2::geom_hline(yintercept = -log10(.data$pval), linetype = 2) +
    ggplot2::geom_vline(xintercept = c(-.data$log2FC, .data$log2FC), linetype = 2) +

    ggplot2::xlab("log2 fold change") + ggplot2::ylab("-log10 p-value") +
    ggplot2::scale_colour_manual(values=c("darkgreen", "gray", "red") ) +

    ggrepel::geom_label_repel(data=Genes, aes(label=.data$GeneSymbol),
                     # Add extra padding around each text label.
                     box.padding = 0.5,
                     # Add extra padding around each data point.
                     point.padding = 0.5,
                     # Color of the line segments.
                     segment.color = 'black',
                     # Width of the line segments.
                     segment.size = 0.5,
                     # Draw an arrow from the label to the data point.
                     #arrow = arrow(length = unit(0.01, 'npc')),
                     # Strength of the repulsion force.
                     force = 1,
                     # Maximum iterations of the naive repulsion algorithm O(n^2).
                     max.iter = 3e3)

  p
}


#' Annotate DE Table
#' @description Useful function tofilter the excel file or to do a volcano plot coloring the DE and non DE genes
#'              requires a gene expression table with log2FC and padjust column for each gene,
#'              like the one obtained from `DESeq()` call.
#'
#' @param df A Dataframe obtained from a differential expression test like `DESeq2::DESeq()`.
#' @param log2FC A numeric value specifying log2FoldChange threshold
#' @param padjust A numeric value specifying padjusted threshold
#'
#' @return A Dataframe annotated with Upregulated and Downregulated genes. Thresholds are used defined
#' @export
#'
Annot_DE <- function(df, log2FC = 1, padjust = 0.1) {

  df <- data.frame(df)
  df$DEG <- "NotDE"
  df$DEG[which(df$log2FoldChange >= log2FC & df$padj <= padjust)] <- "Upregulated" #Annotate UPregulated genes
  df$DEG[which(df$log2FoldChange <= -log2FC & df$padj <= padjust)] <- "Downregulated" #Annotate DOWNregulated genes

  return(df)
}
