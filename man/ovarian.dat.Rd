\name{ovarian.dat}
\alias{ovarian.dat}
\docType{data}
\title{
Example input data for LumiWCluster.
}
\description{
This is a subset of ovarian methylation data from Houshdaran et al. (2009). 
}
\usage{data(ovarian.dat)}
\format{
  The format is:
List of 4
 $ Probe_ID   : Factor w/ 1505 levels "AATK_E63_R","AATK_P519_R",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ Sample_ID  : chr [1:6] "Cell Line" "Cell Line" "Endometrioid" "Endometrioid" ...
 $ beta       : num [1:1505, 1:42] 9.87e-01 9.70e-01 9.41e-01 1.50e-08 2.10e-08 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:42] "CL_2008.3" "CL_A1847.3" "CL_A2780.3" "CL_CaOV3.3" ...
 $ det.p.value: num [1:1505, 1:42] 3.68e-38 3.68e-38 1.22e-10 1.14e-33 6.34e-18 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:42] "CL_2008.1" "CL_A1847.1" "CL_A2780.1" "CL_CaOV3.1" ...
}
\source{Houshdaran et al. (2009). DNA Methylation Profiles of Ovarian Epithelial Carcinoma Tumors and Cell Lines \emph{PLoS ONE}. \bold{5}2:e9359.}
\references{P.F. Kuan, S. Wang, H. Chu (2010). A statistical framework for Illumina DNA methylation array. \emph{Technical Report}.}
\examples{
data(ovarian.dat)
str(ovarian.dat)
}
\keyword{datasets}
