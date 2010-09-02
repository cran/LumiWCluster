\name{plotBAnorm}
\alias{plotBAnorm}
\title{
Boxplots before and after normalization (internal function)
}
\description{
An internal function to be used by  \code{\link{normalize.GoldenGate}}.
}
\usage{
plotBAnorm(ori.beta, norm.beta, design.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ori.beta}{
original beta values.
}
  \item{norm.beta}{
normalized beta values. 
}
  \item{design.file}{
design file for GoldenGate Cancer Panel I
}
}
\value{
  Returns a figure with 4 panels, i.e., unnormalized and normalized beta against sequence length and GC content. 
}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\author{Pei Fen Kuan}
\seealso{ \code{\link{normalize.GoldenGate}} }

\examples{
data(GoldenGateDesignFile)
data(ovarian.dat)  
norm.dat <- normalize.GoldenGate(Probe_ID = ovarian.dat$Probe_ID, beta = ovarian.dat$beta, 
	det.p.value = ovarian.dat$det.p.value, GoldenGateDesignFile, plot.option = TRUE)
}
