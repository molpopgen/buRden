#include <stat_allstats.hpp>
#include <esm.hpp>

#include <algorithm>

using namespace Rcpp;
using namespace std;

stat_allstats::stat_allstats( const unsigned & nrows,
			      const Rcpp::IntegerVector & ccstatus,
			      const unsigned & esm_k_value,
			      const bool & normalize_calpha,
			      const bool & simplecounts_calpha ) : __chisq(),
								   __calpha(ccstatus,normalize_calpha,simplecounts_calpha),
								   __MB(nrows,count(ccstatus.begin(),ccstatus.end(),1),ccstatus),
								   esmK(esm_k_value)
{
}

void stat_allstats::update()
{
  __chisq.update();
  __calpha.update();
  __MB.update();
}
 
void stat_allstats::operator()(const int & genotype,
			       const int & ccstatus)
{
  __chisq(genotype,ccstatus);
  __calpha(genotype,ccstatus);
  __MB(genotype,ccstatus);
}
 
Rcpp::List stat_allstats::values()
{
  List MBvalues = __MB.values();
  return List::create( Named("esm.stat") = esm( as<NumericVector>(__chisq.values()["values"]), esmK ),
		       Named("esm.K") = esmK,
		       Named("calpha.stat") = as<double>(__calpha.values()["statistic"]),
		       Named("MB.general.stat") = as<double>(MBvalues["general"]),
		       Named("MB.recessive.stat") = as<double>(MBvalues["recessive"]),
		       Named("MB.dominant.stat") = as<double>(MBvalues["dominant"])
		       );
}
