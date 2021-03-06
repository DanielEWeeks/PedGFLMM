% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega2PedGFLMM.R
\name{init_PedGFLMM}
\alias{init_PedGFLMM}
\title{load Mega2 SQLite database and perform initialization for PedGFLMM usage}
\usage{
init_PedGFLMM(db = NULL, verbose = FALSE, traitname = "default")
}
\arguments{
\item{db}{specifies the path of a \bold{Mega2} SQLite database containing study data.}

\item{verbose}{TRUE indicates that diagnostic printouts should be enabled.
This value is saved in the returned environment.}

\item{traitname}{Name of the affection status trait to use to set the case/control values; by default, "default"}
}
\value{
"environment" containing data frames from an SQLite database and some computed values.
}
\description{
This populates the \bold{R} data frames from the specified \bold{Mega2} SQLite database.
}
\note{
\emph{init_PedGFLMM} sets up the schaidPed and pedPer data frames that are used later in the \emph{DOPedGFLMM} calculation.
 In addition, it initializes a matrix to aid
  in translating a genotype allele matrix to a genotype count matrix.

 It also initializes the results data frame \emph{envir$PedGFLMM_results} to zero rows.
}
\examples{
db = system.file("exdata", "seqsimmGFLMM.db", package="PedGFLMM")
ENV = init_PedGFLMM(db, traitname = "default")
ls(ENV)

}
\seealso{
\code{\link{DOPedGFLMM}}, \code{\link{Mega2PedGFLMM}}
}
