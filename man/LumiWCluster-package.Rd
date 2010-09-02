\name{LumiWCluster-package}
\alias{LumiWCluster-package}
\docType{package}
\title{A weighted model based clustering for Illumina Methylation BeadArray
}
\description{
This package implements a weighted model based clustering by incorporating a quality measure score, e.g., detection p-values for Illumina Methylation BeadArray. If the platform design is GoldenGate Cancer Panel I, a normalization function is also provided.
}
\details{
\tabular{ll}{
Package: \tab LumiWCluster\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-06-06\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
}
\author{
Pei Fen Kuan

Maintainer: Pei Fen Kuan <pfkuan@bios.unc.edu>
}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\keyword{ package }
\seealso{ \code{\link{normalize.GoldenGate}}, \code{\link{LumiWCluster}}}

\examples{
data(GoldenGateDesignFile)
data(ovarian.dat)  
norm.dat <- normalize.GoldenGate(Probe_ID = ovarian.dat$Probe_ID, beta = ovarian.dat$beta, 
	det.p.value = ovarian.dat$det.p.value, GoldenGateDesignFile, plot.option = TRUE)

fit.LumiWCluster <- LumiWCluster(norm.dat$norm.logitbeta, norm.dat$w, 
		norm.dat$g, K=c(2:3), lambda=c(0.1,1), adaptive=TRUE, 
		center=TRUE, normalize=FALSE, epsilon=1e-4, max.iter=2000, trace=FALSE)

}
