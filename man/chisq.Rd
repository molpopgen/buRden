% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{chisq}
\alias{chisq}
\title{Chi-squared statistic for a 2x2 table}
\usage{
chisq(a, b, c, d, yates)
}
\arguments{
\item{a}{An observation}

\item{b}{An observation}

\item{c}{An observation}

\item{d}{An observation}

\item{yates}{Apply continuity correction?}
}
\value{
The equivalent of chisq.test( matrix(c(a,b,c,d),nrow=2,byrow=T),correct=yates )$statistic
}
\description{
Chi-squared statistic for a 2x2 table
}
\details{
Calculated internally on a log10 scale.
}

