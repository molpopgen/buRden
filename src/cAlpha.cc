#include <Rcpp.h>
#include <randWrapper.hpp>
#include <algorithm>
#include <map>

using namespace Rcpp;
using namespace std;

//' The c-alpha statistic
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param normalize Return the statistic divided by the square root of its variance.
//' @param simplecounts See Details.
//' @return The c-alpha test statistic.  If normalize = TRUE, then T/sqrt(Z) is returned, otherwise T is returned.
//' @details  When simplecounts = FALSE, heterozygous and homozygous genotypes are treated as different numbers of observations
//' of the mutation.  In other wordes, simplecounts = FALSE is equivalent to colSums( ccdata[status==1,] ).  When simplecounts=TRUE,
//' all nonzero genotype values are treated as the value 1, equivalent to  apply(data[status==1,], 2, function(x) sum(x>0, na.rm=TRUE)).
//' The latter method is used by the R package AssotesteR.  
//' @references Neale, B. M., Rivas, M. A., Voight, B. F., Altshuler, D., Devlin, B., Orho-Melander, M., et al. (2011). Testing for an Unusual Distribution of Rare Variants. PLoS Genetics, 7(3), e1001322. doi:10.1371/journal.pgen.1001322
//' @examples
//' data(rec.ccdata)
//' status = c( rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases) )
//' #get minor allele freqs in 
//' rec.ccdata.MAFS = colSums( rec.ccdata$genos[which(status==0),] )/(2*rec.ccdata$ncontrols)
//' rec.ccdata.calpha = cAlpha(rec.ccdata$genos[,which(rec.ccdata.MAFS <= 0.05)],status)
// [[Rcpp::export]]
double cAlpha( const IntegerMatrix & data,
	       const IntegerVector & status,
	       const bool & normalize = false,
	       const bool & simplecounts = false)
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
	  n_i += (simplecounts) ? (data(ind,site)>0) : data(ind,site);
	  switch( status[ind] )
	    {
	    case 0:
	      break;
	    case 1:
	      //See note above re: AssotesteR
	      y_i += (simplecounts) ? (data(ind,site)>0) : data(ind,site);
	      break;
	    default:
	      stop("cAlpha: phenotype label other than 0 or 1 encountered");
	    }
	}
      T += ( pow( double(y_i)-double(n_i)*p0, 2.) - double(n_i)*p0*(1.-p0) );
      if (normalize)
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
    }

  if( normalize )
    {
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

  return T;
}

//' Permutation distribution of the c-alpha statistic
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param nperms The number of permutations to perform
//' @param simplecounts See Details.
//' @return The distribution of the test statistic after nperms swapping of case/control labels
//' @references Neale, B. M., Rivas, M. A., Voight, B. F., Altshuler, D., Devlin, B., Orho-Melander, M., et al. (2011). Testing for an Unusual Distribution of Rare Variants. PLoS Genetics, 7(3), e1001322. doi:10.1371/journal.pgen.1001322
//' @details  When simplecounts = FALSE, heterozygous and homozygous genotypes are treated as different numbers of observations
//' of the mutation.  In other wordes, simplecounts = FALSE is equivalent to colSums( ccdata[status==1,] ).  When simplecounts=TRUE,
//' all nonzero genotype values are treated as the value 1, equivalent to  apply(data[status==1,], 2, function(x) sum(x>0, na.rm=TRUE)).
//' The latter method is used by the R package AssotesteR.  
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
			   const unsigned & nperms, const bool & simplecounts = false  )
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
