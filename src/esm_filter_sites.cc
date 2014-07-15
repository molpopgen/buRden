#include <esm_filter_sites.hpp>
#include <Sequence/Correlations.hpp>
#include <algorithm>
#include <numeric>

using namespace std;
using namespace Sequence;

vector<short> filter_sites(const CCblock & ccdata,const double & minfreq,const double & maxfreq, const double & rsq_cutoff)
{
  //step 1, filter on frequency in controls
  vector<short> keep(ccdata.SN+ccdata.SC,1);
  
  for( unsigned site_i = 0 ; site_i < ccdata.geno_matrix.size() - 1 ; ++site_i )
    {
      double maf_i = double(accumulate(ccdata.geno_matrix[site_i].begin(),ccdata.geno_matrix[site_i].begin()+ccdata.ncontrols,0))/(2.*double(ccdata.ncontrols));

      if( maf_i < minfreq || maf_i > maxfreq )
	{
	  keep[site_i]=0;
	}
      if ( keep[site_i] )
	{
	  for( unsigned site_j = site_i+1 ; site_j < ccdata.geno_matrix.size()  ; ++site_j )  
	    {
	      double maf_j = double(accumulate(ccdata.geno_matrix[site_j].begin(),ccdata.geno_matrix[site_j].begin()+ccdata.ncontrols,0))/(2.*double(ccdata.ncontrols));
	      if( maf_j < minfreq || maf_j > maxfreq )
		{
		  keep[site_j]=0;
		}
	      if ( keep[site_j] )
		{
		  double corr = ProductMoment()(ccdata.geno_matrix[site_i].begin(),ccdata.geno_matrix[site_i].end(),
						ccdata.geno_matrix[site_j].begin());
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
