#include <filter_sites.hpp>
#include <ProductMoment.hpp>
#include <algorithm>
#include <numeric>

using namespace std;
using namespace Rcpp;

//' Apply frequency and LD filters to a genotype matrix
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param minfreq A site with minor allele frequency < minfreq will not be kept.
//' @param maxfreq A site with minor allele frequency > maxfreq will not be kept.
//' @param rsq_cutoff  When comparing two sites, if the genotype correlation coefficient r^2 is >= rsq_cutoff, only the first site will be kept.
//' @return A vector of integers containing the values 0 (not kept) and 1 (kept).  The length of the vector is equal to the number of columns in ccdata.
//' @details Regarding rsq_cutoff, when sites i and j are compared (j > i), site i will be kept and site j will not be kept.
//' @examples
//' data(rec.ccdata)
//' status=c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' keep=filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
// [[Rcpp::export]]
Rcpp::IntegerVector filter_sites(const Rcpp::IntegerMatrix & ccdata,
				 const Rcpp::IntegerVector & ccstatus,
				 const double & minfreq,
				 const double & maxfreq,
				 const double & rsq_cutoff)
{
  if( ccstatus.size() != ccdata.nrow() )
    {
      stop("filter_sites: length(ccstatus) != nrow(ccdata)");
    }
  //step 1, filter on frequency in controls
  Rcpp::IntegerVector keep(ccdata.ncol(),1);

  unsigned ncontrols = count( ccstatus.begin(),ccstatus.end(),0 );
  for( unsigned site_i = 0 ; site_i < ccdata.ncol() - 1 ; ++site_i )
    {
      if ( keep[site_i] )
	{
	  double maf_i = 0;
	  for( unsigned i = 0 ; i < ccdata.nrow() ; ++i )
	    {
	      if( ccstatus[i] == 0 )//is control
		{
		  maf_i += double(ccdata(i,site_i));
		}
	    }
	  maf_i /= (2*double(ncontrols));
	  if( maf_i < minfreq || maf_i > maxfreq )
	    {
	      keep[site_i]=0;
	    }
	  if ( keep[site_i] )
	    {
	      for( unsigned site_j = site_i+1 ; site_j < ccdata.ncol()  ; ++site_j )  
		{
		  if( keep[site_j] )
		    {
		      double maf_j=0;
		      for( unsigned i = 0 ; i < ccdata.nrow() ; ++i )
			{
			  if( ccstatus[i] == 0 )//is control
			    {
			      maf_j += double(ccdata(i,site_j));
			    }
			}
		      maf_j /= (2*double(ncontrols));
		      
		      if( maf_j < minfreq || maf_j > maxfreq )
			{
			  keep[site_j]=0;
			}
		      if ( keep[site_j] )
			{
			  iterator_traits<NumericVector::const_iterator>::value_type corr = 
			    std::pow(ProductMoment( NumericVector( ccdata(_,site_i).begin(),ccdata(_,site_i).end() ),
						    NumericVector( ccdata(_,site_j).begin(),ccdata(_,site_j).end() ) ), 2.);
			  if( corr > rsq_cutoff )
			    {
			      keep[site_j]=0;
			    }
			}
		    }
		}
	    }
	}
    }
  return keep;
}
