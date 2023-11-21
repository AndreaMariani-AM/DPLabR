#' Quantify ChIPseq Scores
#' @description Custom function to quantify ChIPseq scores in multiple .bw files
#'        across multiple genomics regions in BED format
#' @param bw_list A named list of .bw files
#' @param regions A named list of genomic regions in BED file
#' @param operation Operation to compute from \code{\link[megadepth]{get_coverage}} function. Defaul = "mean"
#'
#' @importFrom megadepth get_coverage
#' @importFrom dplyr mutate bind_rows
#' @importFrom purrr map imap
#'
#' @return A tibble with quantifications
#' @export
#'
get_coverage_custom <- function(bw_list, regions, operation) {


  # Check that bw and annot lists are named
  if(is.null(names(bw_list))) { stop("Bigwig list have to be a named list.")}
  if(is.null(names(regions))) { stop("Annotation list have to be a named list.")}

  bw_list %>%
    purrr::map(quantify_regions, regions, operation) %>%
    purrr::imap(~ dplyr::mutate(.x, sample = .y)) %>%
    dplyr::bind_rows()

}

#' Helper Function for get_coverage_custom
#'
#' @param bw A .bw file
#' @param regions A .BED file
#' @param operation Operation for the computation
#'
#' @return A helper function
#' @noRd
quantify_regions <- function(bw, regions, operation = "mean") {

  regions %>%
    purrr::map(~ megadepth::get_coverage(bigwig_file = bw, op = operation, annotation = .x) ) %>%
    purrr::imap(~ dplyr::mutate(.x, group = .y)) %>%
    purrr::map_dfr(as.data.frame)

}
