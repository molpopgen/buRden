#include <Rcpp.h>
#include <stat_allstats.hpp>
#include <stat_calculator.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
List allBurdenStats( const IntegerMatrix & ccdata,
		     const IntegerVector & ccstatus,
		     const unsigned & esm_K,
		     const bool normalize_calpha = false,
		     const bool simplecount_calpha = false )
{
  stat_allstats f(ccdata.nrow(),ccstatus,esm_K,normalize_calpha,simplecount_calpha);
  return stat_calculator(ccdata,ccstatus,f);
}
