#include <stat_base.hpp>

using namespace Rcpp;

List stat_calculator(const IntegerMatrix & data,
		     const IntegerVector & status,
		     stat_base & f)
{
  IntegerMatrix::const_iterator itr = data.begin(),end=data.end();
  const unsigned nr = data.nrow();
  while(itr!=end)
    {
      unsigned i=0;
      while(i<nr)
	{
	  f( *itr++, status[i++] );
	}
      f.update();
    }
  // for( unsigned site = 0 ; site < data.ncol() ; ++site )
  //   {
  //     for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
  // 	{
  // 	  f( data(ind,site), status[ind] );
  // 	}
  //     f.update();
  //   }
  return f.values();
}
