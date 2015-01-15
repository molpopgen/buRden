#' Turn genotypes into normalized genotypes
#' @param genotypes An individuals x mutation matrix, coded 0,1,2 = number of copies of minor alleles
#' @note copied from Golan's original code at https://sites.google.com/site/davidgolanshomepage/software/pcgc
norm.genos <- function(genotypes)
    {
        ps <- apply(genotypes,2,mean)/2
	return(t(apply(genotypes,1, function(x) (x-2*ps)/sqrt(2*ps*(1-ps)))))
    }

#' Equation 5 from www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#' @param genotypes An individuals x mutation matrix, coded 0,1,2 = number of copies of minor alleles
#' @note copied from Golan's original code at https://sites.google.com/site/davidgolanshomepage/software/pcgc
#' @references Golan et al. (2014) Measuring missing heritability: Inferring the contribution of common variants. www.pnas.org/cgi/doi/10.1073/pnas.1419064111
calc.Gij <- function(genotypes)
    {
        ps <- apply(genotypes,2,mean)/2
	good <- ps!=0 & ps!=1
	genotypes <- genotypes[,good]
	m <- ncol(genotypes)
	ps <- ps[good]
	new.genotypes <- norm.genos(genotypes)
	G <- new.genotypes %*% t(new.genotypes) / m 
	return(G)
    }

#' Write a case/control block to files that will be used
#' as the input to the pcgc software from https://sites.google.com/site/davidgolanshomepage/software/pcgc,
#' which implements the methods in www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#' @param genotypes A matrix of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal nrow(x)
#' @param filenamebase  The prefix to use for output file names.  See note below
#' @references Golan et al. (2014) Measuring missing heritability: Inferring the contribution of common variants. www.pnas.org/cgi/doi/10.1073/pnas.1419064111
#'
#' @examples
#' \dontrun{
#' data(rec.ccdata)
#' rec.ccdata.status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
#' cc2pcgc(rec.ccdata$genos,rec.ccdata.status,"rec.ccdata.pcgc")
#' }
cc2pcgc <- function( genotypes, status, filenamebase )
    {
        #write the raw genos file
        write.table(genotypes,file=filenamebase,col.names=FALSE,row.names=FALSE)
        #write the phenos file
        write.table( cbind(1:nrow(genotypes),1:nrow(genotypes),status),
                    row.names=FALSE,col.names=FALSE,
                    file = paste(filenamebase,".phen",sep="") )

        #Now, the normalized genotype matrix
        genotypes.Gij = calc.Gij(genotypes)
        MM = matrix(10000,choose(nrow(genotypes),2) + nrow(genotypes),4)
        C=1
        for( i in 1:nrow(genotypes) )
            {
                for(j in 1:i)
                    {
                        MM[C,1] <- i
                        MM[C,2] <- j
                        MM[C,4] <- genotypes.Gij[i,j]
                        C <- C + 1
                    }
            }
        ofile = gzfile( paste(filenamebase,".grm.gz",sep=""),"w" )
        write.table( MM, ofile, col.names=F,row.names=F,sep="\t" )
        flush(ofile)
        close(ofile)
        #Finally, the normalized genotype id file
        write.table( cbind(1:nrow(genotypes),1:nrow(genotypes)),
                    file = paste(filenamebase,".grm.id",sep=""),
                    row.names=FALSE,col.names=FALSE)
    }
