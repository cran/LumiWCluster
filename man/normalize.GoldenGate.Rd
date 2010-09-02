\name{normalize.GoldenGate}
\alias{normalize.GoldenGate}
\title{
Normalization for Illumina GoldenGate Cancer Panel I
}
\description{
This function normalizes the data from Illumina GoldenGate Methylation platform by adjusting for sequence length and GC content.
}
\usage{
normalize.GoldenGate(Probe_ID, beta, det.p.value, design.file, plot.option = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Probe_ID}{
a vector of Probe ID. Each entry is the name of the CpG site.
}
  \item{beta}{
a matrix of original beta values. Takes value between 0 and 1. Each row is a CpG, each column is a sample.
}
  \item{det.p.value}{
a matrix of detection p-values from BeadStudio. Each row is a CpG, each column is a sample.
}
  \item{design.file}{
Design file for GoldenGate platform, which is included in this package. See example to load the design file.
}
  \item{plot.option}{
a logical value. If TRUE, boxplots before and after normalization will be generated.
}
}
\value{
   return a list of values.
   \item{norm.beta}{normalized beta values.}
   \item{norm.logitbeta}{normalized logit transformed beta values. To be used in \code{\link{LumiWCluster}}}
   \item{w}{a vector of weights for each sample computed from detection p-values. To be used in \code{\link{LumiWCluster}}}
   \item{g}{a vector of weights for each CpG computed from detection p-values. To be used in \code{\link{LumiWCluster}}} 
}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\author{Pei Fen Kuan}
\seealso{ \code{\link{LumiWCluster}} }

\examples{
data(GoldenGateDesignFile)
data(ovarian.dat)  
norm.dat <- normalize.GoldenGate(Probe_ID = ovarian.dat$Probe_ID, beta = ovarian.dat$beta, 
	det.p.value = ovarian.dat$det.p.value, GoldenGateDesignFile, plot.option = TRUE)

}
