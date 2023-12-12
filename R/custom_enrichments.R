#-------------------------------------------------------------------------------------------------------------------------#
# GENE ENRICHMENT ANALYSIS - Run the previous functions
#-------------------------------------------------------------------------------------------------------------------------#
# The input for this function must be a dataframe with 2 columns:
# One called ENTREZID and another called SYMBOL. Reference can be either "human" or "mouse"
#' Run Enrichments
#' @description This function execute different enrichment analysis, namely GO, PA and 2 msig. Takes in a 2 column dataframe
#'              one called `ENTREZID` and the other `SYBOL`. These two columns represent the mapping between gene names and ID.
#'              This DF can be obtained in different ways, for example with `bitr`
#'
#' @param genes_df A dataframe with two columns `ENTREZID` and `SYMBOL`
#' @param reference A character specifying with reference needs to be used. For now it supports "human" and "mouse"
#'
#' @import clusterProfiler
#' @importFrom rlang .data
#' @importFrom dplyr select
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import ReactomePA
#' @import msigdbr
#'
#' @return A data list of dataframe, one per enrichment,
#' @export
#'
enrichments <- function(genes_df, reference) {

  if(reference == "human"){
    species_set <- "Homo sapiens"
    kegg.genome_set <- "hsa"
    pa.genome_set <-  "human"
    db_set <- "org.Hs.eg.db"
  }
  else{
    species_set <- "Mus musculus"
    kegg.genome_set <- "mmu"
    pa.genome_set <-  "mouse"
    db_set <- "org.Mm.eg.db"
  }

  # Get gene sets to perform enrichment analysis
  #Hallmark
  m_df.h = msigdbr::msigdbr(species = species_set, category = "H") %>%
    dplyr::select(.data$gs_name, .data$human_gene_symbol) %>%
    as.data.frame()

  #C2
  m_df.c2 = msigdbr::msigdbr(species = species_set, category = "C2") %>%
    dplyr::select(.data$gs_name, .data$human_gene_symbol) %>%
    as.data.frame()

  # Set mouse references for the databases
  kegg.genome <- kegg.genome_set
  pa.genome   <- pa.genome_set
  db          <- db_set

  #Run custom functions for enrichments
  go      <- goEnrichment(genes_df, ont = "BP", db = db_set)
  # kegg    <- KEGGenrichment(genes_df, org = kegg.genome, db = db)
  pa      <- PAenrichment(genes_df, org = pa.genome_set)
  msig_db_hallmark <- msig_db_enrichment(genes_df, pathways.gmt = m_df.h)
  msig_db_c2all    <- msig_db_enrichment(genes_df, pathways.gmt = m_df.c2)

  #save each output to a list
  list_of_datasets <- list( "GO"        = go,
                            "Reactome"  = pa,
                            "msig_db_hallmark" = msig_db_hallmark,
                            "msig_db_c2all"    = msig_db_c2all
  )

  return(list_of_datasets)
}


#' Custom GSEA Plot
#'
#' @description This function creates a custom GSEA plot in the form a volcano where all the genesets tested are reported, provided that they
#'              pass a `NES` threshold
#'
#' @param df A dataframe containing the results of differential expression analysis. It assumes the structure we give to those dataframes in our pipelines
#'           and therefore should have column `Geneid` at index 1 and `log2FoldChange` at index 3.
#' @param specie A character. What specie should be used. Default is "human". Anything different than that is evaluated as "mouse"
#' @param query A character. What pathway should be tested. By default every gene set that partially matches that pathway. (e.g "WNT", "APOPTOSIS")
#' @param NES_cutoff A numeric. What Normalized Enrichment Score to use as a cutoff (bidirectional). It doesn't assume any significance at this point.
#'
#' @import magrittr
#' @import dplyr
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import ggplot2
#' @importFrom forcats fct_reorder
#' @importFrom rlang .data
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_detect str_replace_all str_wrap
#' @importFrom tibble remove_rownames deframe
#' @importFrom AnnotationDbi select
#' @importFrom purrr flatten
#' @importFrom glue glue
#' @importFrom msigdbr msigdbr
#' @importFrom fgsea fgsea
#' @importFrom stats na.omit
#'
#' @return A list with 4 elements for convenience.
#' @export
#'
plot_fgsea <- function(df, specie = "human", query, NES_cutoff = 1){

  #Evaluating specie to use
  if (specie == "human"){
    specie <- "Homo sapiens"
    db <- org.Hs.eg.db
  }
  else {
    specie <-  "Mus musculus"
    db <- org.Mm.eg.db
  }


  # Rank genes by descending log2FC and select Geneid
  ranked_genes <- df %>%
    dplyr::arrange(dplyr::desc(.data$log2FoldChange)) %>%
    dplyr::select(1,3)

  #Map Symbols to ENTREZID
  symbol2entrez <- AnnotationDbi::select(db,
                                         key=ranked_genes$Geneid,
                                         columns="ENTREZID",
                                         keytype="SYMBOL") %>%
    dplyr::as_tibble()

  #Joining them together
  gene_list <- dplyr::inner_join(ranked_genes, symbol2entrez, by = c("Geneid"="SYMBOL")) %>%
    dplyr::select(3,2) %>%
    na.omit() %>%
    tibble::remove_rownames()

  #Fix ranked ENTREZ IDs
  ranked_entrez <- tibble::deframe(gene_list)


  # Getting pathways
  msigdbr <- msigdbr::msigdbr(species = specie) %>%
    dplyr::select(.data$gs_cat, .data$gs_subcat, .data$gs_name, .data$gene_symbol, .data$entrez_gene, .data$gs_description)

  # Add to query character
  to_query <- c(query)

  #Retrieve pathways to test
  pathways_to_test <- msigdbr %>%
    dplyr::filter(stringr::str_detect(.data$gs_name, paste(to_query, collapse = "|"))) %>%
    #filter(gs_cat == "H" | gs_cat =="C2" ) %>%
    dplyr::select(.data$gs_name, .data$entrez_gene) %>%
    tidyr::pivot_wider(names_from = .data$gs_name, values_from = .data$entrez_gene, values_fn = list) %>%
    purrr::flatten()


  # Gsea + tidy output
  gsea <- fgsea::fgsea(pathways = pathways_to_test,
                stats = ranked_entrez,
                minSize = 15,
                maxSize = 2000) %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(dplyr::desc(.data$NES))

  #Plot GSEA code chunk
  plot_gsea <- gsea %>%
    dplyr::select(-c(.data$leadingEdge, .data$ES)) %>%
    dplyr::filter(.data$NES < -NES_cutoff | .data$NES > NES_cutoff) %>%
    dplyr::mutate(pathway = stringr::str_replace_all(.data$pathway, pattern = "_", replacement = " ")) %>%
    dplyr::mutate(pathway = stringr::str_wrap(.data$pathway,width=50, whitespace_only = TRUE)) %>%
    dplyr::arrange(.data$NES) %>%
    dplyr::mutate(pathway = forcats::fct_reorder(.data$pathway, .data$NES)) %>%
    ggplot2::ggplot(aes(x = .data$NES, y = -log10(.data$padj), col = .data$padj<0.05)) +
    ggplot2::geom_point(size=6, alpha=0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    ggplot2::geom_vline(xintercept = c(-NES_cutoff,NES_cutoff), linetype="dashed") +
    ggplot2::scale_color_manual(values=c("darkred", "darkblue")) +
    ggplot2::labs(title = glue::glue("Signatures Associated with {query} (NES > {NES_cutoff} or < -{NES_cutoff})"))

  return(list(pathways_to_test = pathways_to_test,
              ranked_entrez = ranked_entrez,
              plot_gsea = plot_gsea,
              gsea = gsea))

}

#' Custom Go Function
#'
#' @description This is the custom function run by the previous analysis. It performs GO analysis correcting with `BH`
#'
#' @param df A 2 column dataframe, one called `ENTREZID` and the other one `SYMBOL` representing genes to be tested.
#' @param ont The ontology to use. Default is `BP`
#' @param db A character specifying with specie database to use. Default is Mouse but will be override by an argument in the
#'            main function \code{\link[DPLabR]{enrichments}}
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom dplyr filter
#'
#' @return A dataframe with ontologies filterd by p.adjusted < 0.05
#' @export
#'
goEnrichment <- function(df, ont = "BP", db = org.Mm.eg.db) {

  ego <- clusterProfiler::enrichGO(gene          = df$ENTREZID,
                                   OrgDb         = db,
                                   keyType       = 'ENTREZID',
                                   ont           = ont,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.1,
                                   readable = TRUE)

  # Return an empty dataframe in case there's not result
  if(is.null(ego)) {
    return(data.frame())
  }
  else{
    return(ego@result %>% dplyr::filter(.data$p.adjust < 0.05) %>% add_row())
  }
}


#' Custom KEGG function
#'
#' @description This is the custom function run by the previous analysis. It performs KEGG analysis correcting with `BH`
#'
#' @param df A 2 column dataframe, one called `ENTREZID` and the other one `SYMBOL` representing genes to be tested.
#' @param org A character specifying which specie to use. Default "mmu" but is overridden in the main function.
#' @param db A character specifying with specie database to use. Default is Mouse but will be override by an argument in the
#'            main function \code{\link[DPLabR]{enrichments}}
#'
#' @importFrom clusterProfiler enrichKEGG setReadable
#' @importFrom rlang .data
#' @importFrom dplyr filter
#'
#' @return A dataframe with KEGG pathways filterd by p.adjusted < 0.05
#' @export
#'
KEGGenrichment <- function(df, org = "mmu", db = org.Mm.eg.db) {

  ekgg <- clusterProfiler::enrichKEGG(gene          = df$ENTREZID,
                                      organism      = org,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.1)

  # If the output is empty don't try to do setReadable to avoid an error.
  if(is.null(ekgg)) {
    return(data.frame())
  }
  else{
    ekgg <- clusterProfiler::setReadable(ekgg, OrgDb = db, keyType="ENTREZID")
    return(ekgg@result %>% dplyr::filter(.data$p.adjust < 0.05) %>% add_row())
  }
}


#' Custom Pathway function
#'
#' @description This is the custom function run by the previous analysis. It performs Pathway analysis.
#'
#' @param df A 2 column dataframe, one called `ENTREZID` and the other one `SYMBOL` representing genes to be tested.
#' @param org A character specifying with specie organism to use. Default is Mouse but will be override by an argument in the
#'            main function \code{\link[DPLabR]{enrichments}}
#'
#' @importFrom ReactomePA enrichPathway
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @return A dataframe with Pathways analysis filtered by p.adjusted < 0.5
#' @export
#'
PAenrichment <- function(df, org = "mouse") {
  ePA <- ReactomePA::enrichPathway(gene         = df$ENTREZID,
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.1,
                                   organism     = org,
                                   readable     = TRUE)

  # Return an empty dataframe in case there's not result
  if(is.null(ePA)){
    return(data.frame())
  }
  else{
    return(ePA@result %>% dplyr::filter(.data$p.adjust < 0.05) %>% add_row())
  }
}


#' Custom Gene Set Enrichments
#'
#' @description This is the custom function run by the previous analysis. It performs gene set enrichment analysis.
#'
#' @param df A 2 column dataframe, one called `ENTREZID` and the other one `SYMBOL` representing genes to be tested.
#' @param pathways.gmt A dataframe containing pathways/gene sets to be tested
#'
#' @importFrom clusterProfiler enricher
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#'
#' @return A dataframe with gene set analysis filtered by p.adjusted < 0.5
#' @export
#'
msig_db_enrichment <- function(df, pathways.gmt) {

  # preprocessing symbol names
  df$SYMBOL <- toupper(df$SYMBOL)

  GSEA_hyp  <- clusterProfiler::enricher(gene        = df$SYMBOL,
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.1,
                                         TERM2GENE    = pathways.gmt)

  if(is.null(GSEA_hyp)){
    return(data.frame())
  }
  else{
    return(GSEA_hyp@result %>% dplyr::filter(.data$p.adjust < 0.05) %>% add_row())
  }
}


#' Fraction to Number
#'
#' @description util function to transform a fraction to a number from the output of GO enrichment
#' @param x fraction from GO enrichment to transform to a float
#'
#' @return a float
#' @export
#'

fract_to_dec <- function(x) {
  dec <- eval(parse(text=x))
  return(dec)
}
