#' Function to parse aligment log file
#' @description This is a function to parse log information from alignment files. For now it retrieves only
#'  discarded multimapped reads. Needs to be expanded in functionalities.
#'
#' @param x A named List of character vectors with one element for each line
#'
#' @return a Named list with one numerical value per element of the list
#' @export
#'
get_log_multimapped <- function(x){

  tmp <- list()
  for (i in names(x)) {
    tmp[[i]] <- x[[i]][4]
    tmp[[i]] <- strsplit(tmp[[i]], split="[\\(\\)]")[[1]][2]
    tmp[[i]] <- as.numeric(gsub("\\%", "", tmp[[i]]))
  }
  return(tmp)
}
