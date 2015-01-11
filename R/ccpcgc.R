#' Turn genotypes into normalized genotypes
#' @note copied from Golan's original code at https://sites.google.com/site/davidgolanshomepage/software/pcgc,
norm.genos <- function(genotypes)
    {
        ps <- apply(genotypes,2,mean)/2
	return(t(apply(genotypes,1, function(x) (x-2*ps)/sqrt(2*ps*(1-ps)))))
    }

#' Equation 5 from www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#' @param genotypes A matrix of genotypes
#' @note copied from Golan's original code at https://sites.google.com/site/davidgolanshomepage/software/pcgc,
#' @references Golan et al. (2014) Measuring missing heritability: Inferring the contribution of common variants. www.pnas.org/cgi/doi/10.1073/pnas.1419064111
calc.Gij <- function(genotypes)
    {
        ps <- apply(genotypes,2,mean)/2
	good <- ps!=0 & ps!=1
	genotypes <- genotypes[,good]
	m <- ncol(genotypes)
	ps <- ps[good]
	new.genotypes <- normalize.genotypes(genotypes)
	G <- new.genotypes %*% t(new.genotypes) / m 
	return(G)
    }

#' Write a case/control block to files that will be used
#' as the input to the pcgc software from https://sites.google.com/site/davidgolanshomepage/software/pcgc,
#' which implements the methods in www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#' @param genotypes A matrix of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal nrow(x)
#' @param filename base  The prefix to use for output file names
#' @references Golan et al. (2014) Measuring missing heritability: Inferring the contribution of common variants. www.pnas.org/cgi/doi/10.1073/pnas.1419064111
cc2pcgc <- function( genotypes, status, filenamebase )
    {
    }
