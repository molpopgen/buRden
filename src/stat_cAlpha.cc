#include <stat_cAlpha.hpp>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace std;

stat_cAlpha::stat_cAlpha(const Rcpp::IntegerVector & status,
			 const bool & normalize,
			 const bool & simplecounts) : stat_base(),T(0.), p0( double(count( status.begin(), status.end(), 1 ))/double( status.size() ) ),
						      n_i(0),y_i(0),norm(normalize),simple(simplecounts),ns(map<unsigned,unsigned>())
{

}


void stat_cAlpha::update()
{
  T += ( pow( double(y_i)-double(n_i)*p0, 2.) - double(n_i)*p0*(1.-p0) );
  if (norm)
    {
      map<unsigned,unsigned>::iterator itr = ns.find(n_i);
      if( itr == ns.end() )
	{
	  ns[n_i]=1;
	}
      else
	{
	  itr->second++;
	}
    }
  n_i=y_i=0;
}

void stat_cAlpha::operator()(const int & genotype,
			     const int & ccstatus)
{
  /*
    NOTE: AssotesteR does things differently.
    That package would do:
    n_i += (data(ind,site)>0)?1:0;
    
    I argue that AssotesteR is wrong, as het- vs hom-
    genotype is probably kinda important.  For example,
    if 10 aa genotypes are in cases vs. 10Aa genotypes in 
    controls, there are 2x as many observations of the mutation
    in cases, which the above calculation misses.
  */
  int val = (simple) ? (genotype>0) : genotype; 
  n_i += val;
  switch( ccstatus )
    {
    case 0:
      break;
    case 1:
      //See not above re: AssotesteR
      y_i += val;
      break;
    default:
      stop("cAlpha: phenotype label other than 0 or 1 encountered");
    }
}

double stat_cAlpha::Z() const
{
  double Z = 0.;
  if( norm )
    {
      for( map<unsigned,unsigned>::const_iterator itr = ns.begin() ; itr != ns.end() ; ++itr )
	{
	  double n = double(itr->first),m_of_n=double(itr->second);
	  double inner=0.;
	  double np01mp0 = n*p0*(1.-p0);
	  double np0 = n*p0;
	  for( unsigned u = 0 ; u <= n ; ++u )
	    {
	      inner += R::dbinom(u,n,p0,0)*pow((pow(u-np0,2.) - np01mp0),2.);
	    }
	  Z += m_of_n*inner;
	}
    }
  else
    {
      Z = numeric_limits<double>::quiet_NaN();
    }
  return Z;
}

Rcpp::List stat_cAlpha::values()
{
  double __Z = Z();
  return List::create( Named("statistic") = (norm) ? T/sqrt(__Z) : T
		       );
}
