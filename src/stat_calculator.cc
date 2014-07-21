#include <stat_base.hpp>

using namespace Rcpp;

List stat_calculator(const IntegerMatrix & data,
		     const IntegerVector & status,
		     stat_base & f)
{
  for( unsigned site = 0 ; site < data.ncol() ; ++site )
    {
      for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
	{
	  f( data(ind,site), status[ind] );
	}
      f.update();
    }
  return f.values();
}
		     
