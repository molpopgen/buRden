#include <chisq_per_marker.hpp>

#include <cstdlib>
#include <cmath>
#include <cassert>
#ifndef NDEBUG
#include <algorithm>
#endif
#include <gsl/gsl_cdf.h>

using namespace std;

double get_log10_chisq(const unsigned ctable[4])
//returns -log10(p-value chisq w/Yate's correction with DF = 1)
{
  double a = double(ctable[0]);
  double b = double(ctable[2]);
  double c = double(ctable[1]);
  double d = double(ctable[3]);

  double N = a+b+c+d;
  double rv = log10(N)+2.*log10(fabs(a*d-b*c)-N/2.) - ( log10(a+b)+log10(c+d)+log10(b+d)+log10(a+c) );
  if (! isfinite(rv) )
    {
      //then the chisquared is 0, the p-value is 1, and -log10(1) = 0 
      return 0;
    }
  return( -log10(gsl_cdf_chisq_Q(pow(10,rv),1)) );  
}

vector<double> chisq_per_marker( const CCblock * ccdata, 
				 const vector<short> * keep, 
				 const vector<short>  & ccstatus )
{
  vector<double> csqs;
  assert( count(ccstatus.begin(),ccstatus.end(),0) == ccdata->ncontrols );
  assert( count(ccstatus.begin(),ccstatus.end(),1) == ccdata->ncases );
  vector<short>::const_iterator kbegin = keep->begin();
  for(unsigned site = 0 ; site < ccdata->geno_matrix.size() ; ++site)
    {
      if(*(kbegin+site))
	{
	  unsigned ctable[4]; //minor in controls, minor in cases, major in controls, major in cases
	  ctable[0]=ctable[1]=ctable[2]=ctable[3]=0;
	  unsigned ind = 0;
#ifndef NDEBUG
	  unsigned STATUS_CHECK[2];
	  STATUS_CHECK[0]=0;
	  STATUS_CHECK[1]=0;
#endif
	  for( ; ind < ccdata->ncontrols + ccdata->ncases ; ++ind )
	    {
	      unsigned MINOR = (ccstatus[ind] == 0) ? 0 : 1;
	      unsigned MAJOR = (ccstatus[ind] == 0) ? 2 : 3;
#ifndef NDEBUG
	      ++STATUS_CHECK[ ccstatus[ind] ] ;
#endif
	      switch( ccdata->geno_matrix[site][ind] )
		{
		case 0: //homozygous major
		  ctable[MAJOR]+=2;
		  break;
		case 1: //major/minor het
		  ctable[MINOR]++;
		  ctable[MAJOR]++;
		  break;
		case 2: //homozygous minor
		  ctable[MINOR]+=2;
		  break;
		default:
		  abort();
		}
	    }
	  /*
	  cerr << ind << ' ' << ctable[0]+ctable[2] << ' ' << 2*ccdata->ncontrols << ' '
	       << ctable[1]+ctable[3] << ' ' << 2*ccdata->ncases << '\n';
	  */
	  assert( STATUS_CHECK[0] == ccdata->ncontrols );
	  assert( STATUS_CHECK[1] == ccdata->ncases );
	  assert( ctable[0]+ctable[2] == 2*ccdata->ncontrols );
	  assert( ctable[1]+ctable[3] == 2*ccdata->ncases );
	  csqs.push_back( get_log10_chisq(ctable) );
	}
    }
  return csqs;
}
