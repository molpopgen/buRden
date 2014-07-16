#include <esm.hpp>
#include <algorithm>
#include <functional>
#include <cmath>

using namespace std;
using namespace Rcpp;

//' Association stat from Thornton, Foran, and Long (2013) PLoS Genetics
//' @param scores A vector of single-marker association test scores, on a -log10 scale
//' @param K the number of markers used to calculate ESM_K
//' @return The ESM_K test statistic value
//' @note See http://www.ncbi.nlm.nih.gov/pubmed/23437004 for detail on the test statistic
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' rec.ccdata.chisq = chisq_per_marker(rec.ccdata$genos[,which(keep==1)],status)
//' rec.ccdata.esm = esm( rec.ccdata.chisq, 50 )
// [[Rcpp::export]]
double esm( const NumericVector & scores, const unsigned & K )
{
  NumericVector s(scores);

  sort(s.begin(),s.end(),greater<double>());

  unsigned ntests = scores.size();
  unsigned nm = min(K,unsigned(scores.size()));

  double rv=0;
  for( unsigned i = 0 ; i < nm ; ++i )
    {
      rv += (s[i] - ( -log10( double(i+1)/double(ntests) ) ) );
    }
  return rv;
}
