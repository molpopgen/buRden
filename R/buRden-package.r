#' buRden.
#'
#' @name buRden
#' @docType package
#' @description An Rcpp package for performing "burden" tests of associations due to rare variants by permutation.
#' @author Kevin R. Thornton \email{krthornt@@uci.edu}
#' @references \url{https://github.com/molpopgen/buRden}
#' @useDynLib buRden, .registration=TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp loadModule Module
#' @examples
#' #Quick-start:
#' data(rec.ccdata)
#' #Generate control/case status vector
#' rec.ccdata.status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
#' #Filter sites: 0 <= MAF in cases < 0.05 && r^2 between pairs < 0.8
#' rec.ccdata.keep = filter_sites(rec.ccdata$genos,rec.ccdata.status,0,5e-2,0.8)
#' #Get stat + permutation p-values
#' NPERMS=1e2
#' rec.ccdata.calpha = calpha.p.perm(rec.ccdata$genos[,which(rec.ccdata.keep==1)],rec.ccdata.status,NPERMS )
#' rec.ccdata.MB= MB.p.perm(rec.ccdata$genos[,which(rec.ccdata.keep==1)],rec.ccdata.status,NPERMS)
#' rec.ccdata.esm = esm.p.perm(rec.ccdata$genos[,which(rec.ccdata.keep==1)],rec.ccdata.status,NPERMS,50)
NULL

#' @title Example case/control data
#' @name rec.ccdata
#' @docType data
#' @format A simulated case/control panel with 3,000 controls plus 3,000 cases.  The data are a list composed of the following elements:
#' genos = the genotypes.  rows = individuals, columns = markers.  Values are number of copies of the derived allele.\cr
#' pos = mutation positions.\cr
#' ncontrols,ncases = number of controls and cases, respectively.  The controls come first in genos.\cr
#' burdens = A matrix.  Each row is the number of causative mutations inherited from each parent.\cr
#' phenos = A matrix representing phenotype values for each individual.  The first column is the genetic component.  The second is the random component.\cr
#' neutral = An integer.  The first neutral columns in genos are markers with an effect size of 0.\cr
#' causative = An integer.  The remaining causative columns following the neutral ones have nonzero effect sizes.\cr
#' @details The data were simulated according the the model described in http://www.ncbi.nlm.nih.gov/pubmed/23437004, using the parameters
#' from that paper.  The mean effect size was 0.1.
"rec.ccdata"

#' @title Effect sizes for example data
#' @name rec.esizes
#' @docType data
#' @format The effect size block corresponding to rec.ccdata.  The format is a matrix of 4 columns, and there is one row for each mutation in the entire population. (In other words, nrow(rec.esize) > ncol(rec.ccdata$genos)).  The columns are:\cr
#' mutation positions\cr
#' mutation effect sizes\cr
#' Number of occurrences of each mutation in the population.  The population size was 2e4 diploids, meaning that the mutation frequency is the value in the matrix divided by 4e4\cr
#' The age of the mutation, in generations.
"rec.esizes"
