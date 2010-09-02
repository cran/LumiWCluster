\name{box_plot_bias}
\alias{box_plot_bias}
\title{Boxplots before and after normalization (internal function)
 
}
\description{
An internal function to be used by  \code{\link{normalize.GoldenGate}}.
}
\usage{
box_plot_bias(x, y, xlabName, ylabName, mainName)
}
\arguments{
  \item{x}{
a vector of covariate to plot, e.g., sequence length or GC content.
}
  \item{y}{
a matrix of beta values.
}
  \item{xlabName}{
label on x-axis.
}
  \item{ylabName}{
label on y-axis.
}
  \item{mainName}{
main title.
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