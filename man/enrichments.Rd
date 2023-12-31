% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_enrichments.R
\name{enrichments}
\alias{enrichments}
\title{Run Enrichments}
\usage{
enrichments(genes_df, reference)
}
\arguments{
\item{genes_df}{A dataframe with two columns \code{ENTREZID} and \code{SYMBOL}}

\item{reference}{A character specifying with reference needs to be used. For now it supports "human" and "mouse"}
}
\value{
A data list of dataframe, one per enrichment,
}
\description{
This function execute different enrichment analysis, namely GO, PA and 2 msig. Takes in a 2 column dataframe
one called \code{ENTREZID} and the other \code{SYBOL}. These two columns represent the mapping between gene names and ID.
This DF can be obtained in different ways, for example with \code{bitr}
}
