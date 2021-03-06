% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_TDAmut_object.R
\name{create_TDAmut_object}
\alias{create_TDAmut_object}
\title{Creates object of class TDAmut to use throughout pipeline}
\usage{
create_TDAmut_object(exp_table, mut_table)
}
\arguments{
\item{exp_table}{transcriptomic data from cohort with format samples x genes}

\item{mut_table}{mutation data from cohort as a data.frame with columns Sample, Gene, Mutation, Type}
}
\value{
Creates a TDAmut object populated with expression and mutation data
}
\description{
Creates object of class TDAmut to use throughout pipeline
}
