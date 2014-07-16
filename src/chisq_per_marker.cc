#include <chisq_per_marker.hpp>

#include <cstdlib>
#include <cmath>
#include <Rmath.h>

using namespace std;
using namespace Rcpp;

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
  return ( -log10( R::pchisq( std::pow(10,rv), 1., 0, 0 )) );
}

//' Single-marker association test based on the chi-squared statistic
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return A vector of -log10(p-values) from a chi-squared test with one degree of freedom.  The chisq test is based on a 2x2 table of minor vs major allele counts in cases vs. controls.
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #Note that the result should be very very similar to logistic regression under additive model...
//' rec.ccdata.chisq = chisq_per_marker(rec.ccdata$genos, status)
// [[Rcpp::export]]
NumericVector chisq_per_marker( const IntegerMatrix & ccdata,
				const IntegerVector & ccstatus )
{
  NumericVector csqs;
  for(unsigned site = 0 ; site < ccdata.ncol() ; ++site)
    {
      unsigned ctable[4]; //minor in controls, minor in cases, major in controls, major in cases
      ctable[0]=ctable[1]=ctable[2]=ctable[3]=0;
      unsigned ind = 0;
#ifndef NDEBUG
      unsigned STATUS_CHECK[2];
      STATUS_CHECK[0]=0;
      STATUS_CHECK[1]=0;
#endif
      for( ; ind < ccdata.nrow() ; ++ind )
	{
	  unsigned MINOR = (ccstatus[ind] == 0) ? 0 : 1;
	  unsigned MAJOR = (ccstatus[ind] == 0) ? 2 : 3;
	  switch( ccdata(ind,site) )
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
	      stop("chisq_per_marker error: genotype value other than 0, 1, or 2 was encountered!\n");
	    }
	}
      csqs.push_back( get_log10_chisq(ctable) );
    }
  return csqs;
}
