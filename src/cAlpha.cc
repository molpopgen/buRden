#include <Rcpp.h>
#include <randWrapper.hpp>
#include <algorithm>
#include <map>

using namespace Rcpp;
using namespace std;

//' The c-alpha statistic
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return The c-alpha test statistic
//' @examples
//' data(rec.ccdata)
//' status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases) )
//' #get minor allele freqs in 
//' rec.ccdata.MAFS = colSums( rec.ccdata$genos[which(status==0),] )/(2*rec.ccdata$ncontrols)
//' rec.ccdata.calpha = cAlpha(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status)
// [[Rcpp::export]]
double cAlpha( const IntegerMatrix & data,
	       const IntegerVector & status )
{
  /*
    p0 = proportion of cases
    
    For a mutation that occurs n_i total times in the sample,
    we expect there to be y_i = po*n_i in cases and (1-p0)*n_i copies in controls
  */
  unsigned ncontrols = count( status.begin(), status.end(), 0 ),
    ncases = status.size() - ncontrols;
  double p0 = double(ncases)/double(status.size());
  map<unsigned,unsigned> ns;
  double T = 0.;
  for( unsigned site = 0 ; site < data.ncol() ; ++site )
    {
      unsigned n_i=0,y_i=0;
      for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
	{
	  n_i += data(ind,site);
	  switch( status[ind] )
	    {
	    case 0:
	      break;
	    case 1:
	      y_i += data(ind,site);
	      break;
	    default:
	      stop("cAlpha: phenotype label other than 0 or 1 encountered");
	    }
	}
      T += ( pow( double(y_i)-double(n_i)*p0, 2.) - double(n_i)*p0*(1.-p0) );
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

  double Z = 0.;
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
  return ( T/sqrt(Z) );
}

//' Permutation distribution of the c-alpha statistic
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param nperms The number of permutations to perform
//' @return The distribution of the test statistic after nperms swapping of case/control labels
//' @examples
//' data(rec.ccdata)
//' status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases) )
//' #get minor allele freqs in 
//' rec.ccdata.MAFS = colSums( rec.ccdata$genos[which(status==0),] )/(2*rec.ccdata$ncontrols)
//' rec.ccdata.calpha = cAlpha(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status)
//' rec.ccdata.calpha.permdist = cAlpha_perm(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status,100)
// [[Rcpp::export]]
NumericVector cAlpha_perm( const IntegerMatrix & data,
			   const IntegerVector & status,
			   const unsigned & nperms )
{
  NumericVector rv(nperms);
  IntegerVector cc = clone(status);

  for(unsigned i = 0 ; i < nperms ;++i )
    {
      RNGScope scope;
      random_shuffle(cc.begin(),cc.end(),randWrapper);
      rv[i]=cAlpha(data,cc);
    }
  return rv;
}
