\name{posctrl.dat}
\alias{posctrl.dat}
\docType{data}
\title{
Positive control data 
}
\description{
This is the positive control SssI2 methylation data from Houshdaran et al. (2009). 
}
\usage{data(posctrl.dat)}
\format{
  The format is:
List of 4
 $ Probe_ID   : Factor w/ 1505 levels "AATK_E63_R","AATK_P519_R",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ Sample_ID  : chr "Pos Ctrl"
 $ beta       : num [1:1505] 0.79 0.78 0.5 0.48 0.97 0.95 0.91 0.79 0.8 0.98 ...
 $ det.p.value: num [1:1505] 2.63e-31 3.68e-38 3.93e-07 2.24e-22 1.52e-24 ...
}
\source{Houshdaran et al. (2009). DNA Methylation Profiles of Ovarian Epithelial Carcinoma Tumors and Cell Lines \emph{PLoS ONE}. \bold{5}2:e9359.}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\examples{
data(posctrl.dat)
str(posctrl.dat)
}
\keyword{datasets}

