##' Single-marker association test based on the chi-squared statistic
##' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
##' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
##' @param yates Apply continuity correction?
##' @return A vector of -log10(p-values) from a chi-squared test with one degree of freedom.  The chisq test is based on a 2x2 table of minor vs major allele counts in cases vs. controls.
##' @examples
##' data(rec.ccdata)
##' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
##' #Note that the result should be very very similar to logistic regression under additive model...
##' rec.ccdata.chisq = chisq_per_marker(rec.ccdata$genos, status)
chisq_per_marker <- function(ccdata,ccstatus,yates=TRUE)
    {
        .chisq_per_marker(ccdata,ccstatus,yates)
    }
