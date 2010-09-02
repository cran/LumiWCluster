\name{GoldenGateDesignFile}
\alias{GoldenGateDesignFile}
\docType{data}
\title{
Design file for Illumina GoldenGate Cancer Panal I
}
\description{
Design file for Illumina GoldenGate Cancer Panel I with sequence length and GC content for each probe.
}
\usage{data(GoldenGateDesignFile)}
\format{
  A data frame with 1505 observations on the following 3 variables.
  \describe{
    \item{\code{Probe_ID}}{a character vector of Probe ID. Must match Probe ID in input files.}
    \item{\code{L}}{a numeric vector representing the sequence length for each probe.}
    \item{\code{GC}}{a numeric vector representing the GC content for each probe.}
  }
}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\examples{
data(GoldenGateDesignFile)
str(GoldenGateDesignFile)
}
\keyword{datasets}
