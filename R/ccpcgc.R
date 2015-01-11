norm.genos <- function(x)
    {
    }

#' Equation 5 from www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#' @references Golan et al. (2014) Measuring missing heritability: Inferring the contribution of common variants. www.pnas.org/cgi/doi/10.1073/pnas.1419064111
calc.Gij <- function(x)
    {
    }

#' Write a case/control block to files that will be used
#' as the input to the pcgc software from https://sites.google.com/site/davidgolanshomepage/software/pcgc,
#' which implements the methods in www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @param filename base  The prefix to use for output file names
#' @references Golan et al. (2014) Measuring missing heritability: Inferring the contribution of common variants. www.pnas.org/cgi/doi/10.1073/pnas.1419064111
cc2pcgc <- function( x, status, filenamebase )
    {
    }
