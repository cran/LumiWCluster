\name{inv_logit}
\alias{inv_logit}
\title{
Inverse logit function (internal function)
}
\description{
An internal function to be used by  \code{\link{normalize.GoldenGate}}. Computes exp(mu)/(1+exp(mu))
}
\usage{
inv_logit(mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
a vector or matrix of values to be transformed.
}
}
\value{
Returns exp(mu)/(1+exp(mu))
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
