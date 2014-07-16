#include <esm.hpp>
#include <chisq_per_marker.hpp>
#include <randWrapper.hpp>
#include <functional>
#include <utility>
#include <algorithm>

using namespace Rcpp;
using namespace R;
using namespace std;

//' Obtain permutaion distribution of the ESM_K statistic for case/control data
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param nperms Number of permutations to perform
//' @param k Number of markers to use for ESM_K statistic
//' @return A vector of the permuted test statistic values
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' rec.ccdata.chisq = chisq_per_marker(rec.ccdata$genos[,which(keep==1)],status)
//' rec.ccdata.esm = esm( rec.ccdata.chisq, 50 )
//' rec.ccdata.esm.permdist = esm_perm_binary(rec.ccdata$genos[,which(keep==1)],status,100,50)
// [[Rcpp::export]]
NumericVector esm_perm_binary( const IntegerMatrix & ccdata,
			       const IntegerVector & ccstatus,
			       const unsigned & nperms,
			       const unsigned & k )
{
  NumericVector rv(nperms);
  RNGScope scope;
  IntegerVector status = clone(ccstatus);

  for( unsigned i = 0 ; i < nperms ; ++i )
    {
      random_shuffle(status.begin(),status.end(),randWrapper);
      NumericVector c = chisq_per_marker( ccdata, status );
      rv[i] = esm(c,k);
    }
  return rv;
}
		     
