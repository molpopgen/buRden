#include <Rcpp.h>
#include <algorithm>
#include <map>

using namespace Rcpp;
using namespace std;

//' The c-alpha statistic
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @examples
//' data(rec.ccdata)
//' status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases) )
//' rec.ccdata.calpha = cAlpha(rec.ccdata$genos,status)
// [[Rcpp::export]]
double cAlpha( const IntegerMatrix & data,
	       const IntegerVector & status )
{
  unsigned ncontrols = count( status.begin(), status.end(), 0 ),
    ncases = status.size() - ncontrols;
  double p0 = double(ncases)/double(status.size());
  map<unsigned,unsigned> ns;
  double T = 0.;
  for( unsigned site = 0 ; site < data.ncol() ; ++site )
    {
      unsigned co = 0,ca = 0;
      for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
	{
	  switch( status[ind] )
	    {
	    case 0:
	      co += data(ind,site);
	      break;
	    case 1:
	      ca += data(ind,site);
	      break;
	    default:
	      stop("cAlpha: phenotype label other than 0 or 1 encountered");
	    }
	}
      T += ( pow( double(ca)-double(ca+co)*p0, 2.) - double(ca+co)*p0*(1.-p0) );
      unsigned n = ca + co;
      map<unsigned,unsigned>::iterator itr = ns.find(n);
      if( itr == ns.end() )
	{
	  ns[n]=1;
	}
      else
	{
	  itr->second++;
	}
    }
  double Z = 0.;

  for( map<unsigned,unsigned>::const_iterator itr = ns.begin() ; itr != ns.end() ; ++itr )
    {
      Rcerr << itr->first << '\n';
      double inner=0.;
      for( unsigned u = 0 ; u < itr->first ; ++u )
	{
	  inner += R::dbinom(double(u),double(itr->first),p0,0)*pow((pow(double(u)-double(itr->first)*p0,2.) - double(itr->first)*p0*(1.-p0)),2.);
	}
      Z += double(itr->second)*inner;
    }
  Rcerr << T << ' ' << Z << '\n';
  return ( T/sqrt(Z) );
}
