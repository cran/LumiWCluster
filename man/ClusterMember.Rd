\name{ClusterMember}
\alias{ClusterMember}
\title{
Summary of cluster membership
}
\description{
Prints the sample ID for each cluster.
}
\usage{
ClusterMember(LumiWCluster, Sample_ID)
}
\arguments{
  \item{LumiWCluster}{
output object from \code{\link{LumiWCluster}}.}
  \item{Sample_ID}{
lists of sample IDs
}
}
\value{
Return summary printout of sample ID within each cluster.
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

ClusterMember(fit.LumiWCluster,ovarian.dat$Sample_ID)
}