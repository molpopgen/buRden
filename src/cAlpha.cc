#include <Rcpp.h>
#include <stat_cAlpha.hpp>
#include <stat_calculator.hpp>
#include <randWrapper.hpp>
#include <algorithm>
#include <map>

using namespace Rcpp;
using namespace std;

//' The c-alpha statistic
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param normalize Return the statistic divided by the square root of its variance.
//' @param simplecounts See Details.
//' @return The c-alpha test statistic.  If normalize = TRUE, then T/sqrt(Z) is returned, otherwise T is returned.
//' @details  When simplecounts = FALSE, heterozygous and homozygous genotypes are treated as different numbers of observations
//' of the mutation.  In other wordes, simplecounts = FALSE is equivalent to colSums( ccdata[status==1,] ).  When simplecounts=TRUE,
//' all nonzero genotype values are treated as the value 1, equivalent to  apply(data[status==1,], 2, function(x) sum(x>0, na.rm=TRUE)).
//' The latter method is used by the R package AssotesteR.  
//' @references Neale, B. M., Rivas, M. A., Voight, B. F., Altshuler, D., Devlin, B., Orho-Melander, M., et al. (2011). Testing for an Unusual Distribution of Rare Variants. PLoS Genetics, 7(3), e1001322. doi:10.1371/journal.pgen.1001322
//' @examples
//' data(rec.ccdata)
//' status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases) )
//' #get minor allele freqs in the data
//' rec.ccdata.MAFS = colSums( rec.ccdata$genos[which(status==0),] )/(2*rec.ccdata$ncontrols)
//' rec.ccdata.calpha = cAlpha(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status)
// [[Rcpp::export]]
double cAlpha( const IntegerMatrix & ccdata,
	       const IntegerVector & ccstatus,
	       const bool & normalize = false,
	       const bool & simplecounts = false)
{
  stat_cAlpha f(ccstatus,normalize,simplecounts);
  List rv =  stat_calculator(ccdata,ccstatus,f);
  return as<double>(rv["statistic"] );
}

//' Permutation distribution of the c-alpha statistic
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param nperms The number of permutations to perform
//' @param simplecounts See Details.
//' @return The distribution of the test statistic after nperms swapping of case/control labels
//' @references Neale, B. M., Rivas, M. A., Voight, B. F., Altshuler, D., Devlin, B., Orho-Melander, M., et al. (2011). Testing for an Unusual Distribution of Rare Variants. PLoS Genetics, 7(3), e1001322. doi:10.1371/journal.pgen.1001322
//' @details  When simplecounts = FALSE, heterozygous and homozygous genotypes are treated as different numbers of observations
//' of the mutation.  In other wordes, simplecounts = FALSE is equivalent to colSums( ccdata[status==1,] ).  When simplecounts=TRUE,
//' all nonzero genotype values are treated as the value 1, equivalent to  apply(data[status==1,], 2, function(x) sum(x>0, na.rm=TRUE)).
//' The latter method is used by the R package AssotesteR.  
//' @examples
//' data(rec.ccdata)
//' status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases) )
//' #get minor allele freqs in the data
//' rec.ccdata.MAFS = colSums( rec.ccdata$genos[which(status==0),] )/(2*rec.ccdata$ncontrols)
//' rec.ccdata.calpha = cAlpha(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status)
//' rec.ccdata.calpha.permdist = cAlpha_perm(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status,100)
// [[Rcpp::export]]
NumericVector cAlpha_perm( const IntegerMatrix & ccdata,
			   const IntegerVector & ccstatus,
			   const unsigned & nperms, const bool & simplecounts = false  )
{
  NumericVector rv(nperms);
  IntegerVector cc = clone(ccstatus);

  for(unsigned i = 0 ; i < nperms ;++i )
    {
      RNGScope scope;
      random_shuffle(cc.begin(),cc.end(),randWrapper);
      rv[i]=cAlpha(ccdata,cc);
    }
  return rv;
}
