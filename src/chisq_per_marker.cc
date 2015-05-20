#include <chisq_per_marker.hpp>
#include <stat_chisq.hpp>
#include <stat_calculator.hpp>
#include <cstdlib>
#include <cmath>
#include <Rmath.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export(".chisq_per_marker")]]
NumericVector chisq_per_marker( const IntegerMatrix & ccdata,
				const IntegerVector & ccstatus,
				const bool & yates)
{
  stat_chisq f(yates);
  return ( as<NumericVector>(stat_calculator( ccdata, ccstatus, f )["values"]) );
}
