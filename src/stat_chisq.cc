#include <stat_chisq.hpp>
#include <chisq.hpp>
#include <cmath>

using namespace Rcpp;
using namespace std;

stat_chisq::stat_chisq(const bool & use_yates) : stat_base(),yates(use_yates),csqs( NumericVector() )
{
  ctable[0]=ctable[1]=ctable[2]=ctable[3]=0;
}

double stat_chisq::log10chisq()
{
  // double a = double(ctable[0]);
  // double b = double(ctable[2]);
  // double c = double(ctable[1]);
  // double d = double(ctable[3]);

  //double rv = chisq(a,b,c,d,true);
  double rv = chisq(ctable[0],ctable[2],ctable[1],ctable[3],yates);
  if (! isfinite(rv) )
    {
      //then the chisquared is 0, the p-value is 1, and -log10(1) = 0 
      return 0;
    }
  return ( -log10( R::pchisq( rv, 1., 0, 0 )) );
}

void stat_chisq::update() 
{
  csqs.push_back( this->log10chisq() );
  ctable[0]=ctable[1]=ctable[2]=ctable[3]=0;
}

void stat_chisq::operator()(const int & genotype,
			    const int & ccstatus) 
{
  unsigned MINOR = (ccstatus == 0) ? 0 : 1;
  unsigned MAJOR = (ccstatus == 0) ? 2 : 3;
  switch( genotype )
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

List stat_chisq::values()
{
  return List::create(Named("values") = csqs);
}
