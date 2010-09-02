\name{InformativeCpG}
\alias{InformativeCpG}
\title{
List of informative/non-zero CpG loci
}
\description{
This function outputs the list of CpG loci which are not shrunk to zero.
}
\usage{
InformativeCpG(LumiWCluster, Probe_ID)
}
\arguments{
  \item{LumiWCluster}{
output object from \code{\link{LumiWCluster}}.}
  \item{Probe_ID}{
a vector of probe IDs
}
}
\value{
Return a vector of CpG loci which are not shrunk to zero by \code{\link{LumiWCluster}}. These are potentially informative CpGs.
}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\author{Pei Fen Kuan}
\seealso{ \code{\link{LumiWCluster}} }

\examples{
data(GoldenGateDesignFile)
data(ovarian.dat)  
norm.dat <- normalize.GoldenGate(Probe_ID = ovarian.dat$Probe_ID, beta = ovarian.dat$beta, 
	det.p.value = ovarian.dat$det.p.value, GoldenGateDesignFile, plot.option = FALSE)

fit.LumiWCluster <- LumiWCluster(norm.dat$norm.logitbeta, norm.dat$w, 
		norm.dat$g, K=c(2:3), lambda=c(0.1,1), adaptive=TRUE, 
		center=TRUE, normalize=FALSE, epsilon=1e-4, max.iter=2000, trace=FALSE)

Sig_Probes <- InformativeCpG(fit.LumiWCluster, ovarian.dat$Probe_ID)
str(Sig_Probes)
}
