#include <chisq.hpp>
#include <cmath>
#include <algorithm>
//#include <Rcpp.h>

//using namespace Rcpp;

//' Chi-squared statistic for a 2x2 table
//' @param a An observation
//' @param b An observation
//' @param c An observation
//' @param d An observation
//' @param yates Apply continuity correction?
//' @return The equivalent of chisq.test( matrix(c(a,b,c,d),nrow=2,byrow=T),correct=yates )$statistic
//' @details Calculated internally on a log10 scale.
// [[Rcpp::export]]
double chisq(const unsigned & a,
	     const unsigned & b,
	     const unsigned & c,
	     const unsigned & d,
	     const bool & yates)
{
  double _a=a,_b=b,_c=c,_d=d;
  double __N = double(_a+_b+_c+_d);
  double inner = std::max( 0.,fabs(_a*_d-_b*_c) - ( (yates) ? __N/2. : 0 ) );
  if ( inner == 0. ) 
    {
      return 0.;
    }
  double rv = log10(__N)+2.*log10(inner) - ( log10(_a+_b)+log10(_c+_d)+log10(_b+_d)+log10(_a+_c) );
  return std::pow(10,rv);
}
