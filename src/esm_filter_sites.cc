#include <esm_filter_sites.hpp>
#include <ProductMoment.hpp>
#include <algorithm>
#include <numeric>

using namespace std;
using namespace Rcpp;

Rcpp::IntegerVector filter_sites(const Rcpp::IntegerMatrix & ccdata,
				 const Rcpp::IntegerVector & ccstatus,
				 const double & minfreq,
				 const double & maxfreq,
				 const double & rsq_cutoff)
{
  //step 1, filter on frequency in controls
  Rcpp::IntegerVector keep(ccstatus.size(),1);

  unsigned ncontrols = count( ccstatus.begin(),ccstatus.end(),0 );
  for( unsigned site_i = 0 ; site_i < ccdata.ncol() - 1 ; ++site_i )
    {
      double maf_i = 0;
      for( unsigned i = 0 ; i < ccdata.nrow() ; ++i )
	{
	  if( ccstatus[i] == 0 )//is control
	    {
	      maf_i += double(ccdata(i,site_i));
	    }
	}
      maf_i /= 2*double(ncontrols);


      if( maf_i < minfreq || maf_i > maxfreq )
	{
	  keep[site_i]=0;
	}
      if ( keep[site_i] )
	{
	  for( unsigned site_j = site_i+1 ; site_j < ccdata.ncol()  ; ++site_j )  
	    {
	      double maf_j=0;
	      for( unsigned i = 0 ; i < ccdata.nrow() ; ++i )
		{
		  if( ccstatus[i] == 0 )//is control
		    {
		      maf_j += double(ccdata(i,site_j));
		    }
		}
	      maf_j /= 2*double(ncontrols);

	      if( maf_j < minfreq || maf_j > maxfreq )
		{
		  keep[site_j]=0;
		}
	      if ( keep[site_j] )
		{
		  iterator_traits<NumericVector::const_iterator>::value_type corr = 
		    ProductMoment( NumericVector( ccdata(_,site_i).begin(),ccdata(_,site_i).end() ),
				   NumericVector( ccdata(_,site_j).begin(),ccdata(_,site_j).end() ) );
		  if( corr > rsq_cutoff )
		    {
		      keep[site_j]=0;
		    }
		}
	    }
	}
    }
  return keep;
}
