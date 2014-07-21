#include <Rcpp.h>
#include <stat_allstats.hpp>
#include <stat_calculator.hpp>
#include <randWrapper.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace std;

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

// [[Rcpp::export]]
List allBurdenStatsPerm( const IntegerMatrix & ccdata,
			 const IntegerVector & ccstatus,
			 const unsigned & esm_K,
			 const unsigned & nperms,
			 const bool normalize_calpha = false,
			 const bool simplecount_calpha = false )
{
  RNGScope scope;
  IntegerVector status = clone(ccstatus);
  //store permutation distributions
  NumericVector esm_p(nperms),
    calpha_p(nperms),
    MBg_p(nperms),
    MBr_p(nperms),
    MBd_p(nperms);

  for( unsigned i = 0 ; i < nperms; ++i )
    {
      random_shuffle(status.begin(),status.end(),randWrapper); 
      stat_allstats f(ccdata.nrow(),status,esm_K,normalize_calpha,simplecount_calpha);
      List perm_vals = stat_calculator(ccdata,status,f);
      esm_p[i] = as<double>( perm_vals["esm.stat"] );
      calpha_p[i] = as<double>( perm_vals["calpha.stat"] );
      MBg_p[i] = as<double>( perm_vals["MB.general.stat"] );
      MBr_p[i] = as<double>( perm_vals["MB.recessive.stat"] );
      MBd_p[i] = as<double>( perm_vals["MB.dominant.stat"] );
    }
  return List::create( Named("esm.permdist") = esm_p,
		       Named("calpha.permdist") = calpha_p,
		       Named("MB.general.permdist") = MBg_p,
		       Named("MB.recessive.permdist") = MBr_p,
		       Named("MB.dominant.permdist") = MBd_p 
		       );
}
