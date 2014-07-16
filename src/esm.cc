#include <esm.hpp>
#include <algorithm>
#include <functional>
#include <cmath>

using namespace std;
using namespace Rcpp;

double esm( const NumericVector & scores, const unsigned & K )
{
  //vector<double> s(scores);
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
