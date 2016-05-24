\name{km2gcndata}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Gene expression and WGCNA network for 1K genes
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Microarray gene expression and WGCNA network for 1K genes from UKBEC White Matter's microarray.
}
\usage{data("km2gcndata")}
\format{
  A list with two objects.
  \describe{
    \item{\code{kmgcndat$net}}{a WGCNA network}
    \item{\code{kmgcndat$expr}}{a gene expression profiling data frame with genes in columns and samples in rows}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(km2gcndata)
net = applykM2WGCNA(net.label="dummy",
  net.file=km2gcndata$net,
  expr.data=km2gcndata$expr,
  job.path="~/tmp/",
  meg=0)
  str(net)
}
\keyword{datasets}
