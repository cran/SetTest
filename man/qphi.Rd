% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qphi.R
\name{qphi}
\alias{qphi}
\title{Quantile of phi-divergence statitic under the null hypothesis.}
\usage{
qphi(p, M, k0, k1, s = 2, onesided = FALSE)
}
\arguments{
\item{p}{-  a scalar left probability that defines the quantile.}

\item{M}{- correlation matrix of input statistics (of the input p-values).}

\item{k0}{- search range starts from the k0th smallest p-value.}

\item{k1}{- search range ends at the k1th smallest p-value.}

\item{s}{- the phi-divergence test parameter.}

\item{onesided}{- TRUE if the input p-values are one-sided.}
}
\value{
Quantile of phi-divergence statistics.
}
\description{
Quantile of phi-divergence statitic under the null hypothesis.
}
\examples{
## The 0.05 critical value of HC statistic when n = 10:
qphi(p=.95, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE)
}
\references{
1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.

2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
}
\seealso{
\code{\link{stat.phi}} for the definition of the statistic.
}
