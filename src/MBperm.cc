#include <Rcpp.h>
#include <randWrapper.hpp>
#include <MBstat.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace std;

//' Get permutation distribution of Madsen-Browning test statistics
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.                                                                                    
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case. 
//' @param nperms The number of permutations to perform
//' @return A data frame of permuted statistics
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' mbstats = MBstat( rec.ccdata$genos[,which(keep==1)], status )
//' mbstats.perm = MB_perm( rec.ccdata$genos[,which(keep==1)], status, 100 )
// [[Rcpp::export]]
DataFrame MB_perm( const IntegerMatrix & ccdata,
		   const IntegerVector & ccstatus,
		   const unsigned & nperms )
{
  NumericVector g(nperms),d(nperms),r(nperms);
  RNGScope scope;
  IntegerVector status = clone(ccstatus);
 
  for( unsigned i = 0 ; i < nperms ; ++i )
    {
      random_shuffle(status.begin(),status.end(),randWrapper);
      List mbstats = MBstat( ccdata, status );
      g[i] = as<double>(mbstats["general"]);
      d[i] = as<double>(mbstats["dominant"]);
      r[i] = as<double>(mbstats["recessive"]);
    }

  return DataFrame::create( Named("general") = g,
			    Named("recessive") = r,
			    Named("dominant") = d );
}
