% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpecFunctions.R
\name{create_insilico_mixture_template}
\alias{create_insilico_mixture_template}
\title{create_insilico_mixture_template}
\usage{
create_insilico_mixture_template(mono.info, mz.tol = 0.5)
}
\arguments{
\item{mono.info}{An output of \code{\link{gather_summary}}.}

\item{mz.tol}{A m/z tolerance in Da. (Default: 0.5)}
}
\value{
A data frame which contains simulated m/z, log intensity, and normalized intensity values of peaks.
}
\description{
This function generates an intial template for simulated mass spectra.
}
\examples{
spectra.processed.A <- process_monospectra(
   file=system.file("extdata", "listA.txt", package="MGMS2"),
   mass.range=c(1000,2200))
spectra.processed.B <- process_monospectra(
   file=system.file("extdata", "listB.txt", package="MGMS2"),
   mass.range=c(1000,2200))
spectra.processed.C <- process_monospectra(
   file=system.file("extdata", "listC.txt", package="MGMS2"),
   mass.range=c(1000,2200))
spectra.mono.summary.A <- summarize_monospectra(
   processed.obj=spectra.processed.A,
   species='A', directory=tempdir())
spectra.mono.summary.B <- summarize_monospectra(
   processed.obj=spectra.processed.B,
   species='B', directory=tempdir())
spectra.mono.summary.C <- summarize_monospectra(
   processed.obj=spectra.processed.C,
   species='C', directory=tempdir())
mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
template <- create_insilico_mixture_template(mono.info)
}
