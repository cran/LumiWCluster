\name{LumiWCluster}
\alias{LumiWCluster}
\title{
Weighted model based clustering 
}
\description{
This is the main function that implements a weighted model based clustering.
}
\usage{
LumiWCluster(logitbeta, w, g, K = c(2:5), lambda = exp(log(seq(1, 50, length = 20))), 
adaptive = TRUE, center = TRUE, normalize = FALSE, epsilon = 1e-04, max.iter = 2000, 
trace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logitbeta}{
data matrix for clustering taking values on real line. Each row is a CpG, each column is a sample. Remember to logit transform the beta values.
}
  \item{w}{
a vector of weights for each sample. If the data is from Illumina BeadArray platform, the weights can be set as a function of detection p-values.
}
  \item{g}{
a vector of weights for each CpG. If the data is from Illumina BeadArray platform, the weights can be set as a function of detection p-values.
}
  \item{K}{
candidate number of clusters to be considered.
}
  \item{lambda}{
candidate tuning parameters for the penalty function. Larger lambda will shrink more CpGs to zeros.
}
  \item{adaptive}{
a logical value. If TRUE, use adaptive method.
}
  \item{center}{
a logical value. If TRUE, data will be centered.
}
  \item{normalize}{
a logical value. If TRUE, data will be scaled.
}
  \item{epsilon}{
a scalar representing convergence control.
}
  \item{max.iter}{
a scalar representing maximum number of interation. 
}
  \item{trace}{
a logical value. If TRUE, each iteration is printed.
}
}
\value{
return a list of values.
\item{opt.K}{optimal number of clusters.}
\item{opt.lambda}{optimal tuning parameter lambda value.}
\item{clusterID}{cluster membership for each sample.}
\item{BIC}{optimal BIC value.}
\item{BIC_iter}{BIC values for each of candidate number of clusters K and tuning parameter lambda.}
\item{lambda}{candidate tuning parameters for the penalty function.}
\item{K}{candidate number of clusters.}
\item{converge}{if TRUE, the algorithm converges. If FALSE, increase max.iter above.}
\item{pi.vec}{a vector of mixing proportion of the mixture model.}
\item{mu}{a matrix of cluster means for each CpG.}
\item{mu.centering}{a matrix of cluster means for each CpG for the centered data.}
\item{init.mu}{a matrix of initial cluster means for each CpG.}
\item{sigma2.vec}{a matrix of cluster variance for each CpG.}
\item{iter}{number of iterations}

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
}
