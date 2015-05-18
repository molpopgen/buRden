#include<MBstat.hpp>
#include <stat_MadsenBrowning.hpp>
#include <stat_calculator.hpp>
#include <randWrapper.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

using namespace Rcpp;
using namespace std;

//' Calculate Madsen-Browning weights.
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return An array of weights, one for each column in data.
//' @details Calculation is done under the "general genetic model" defined in Madsen and Browning.
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
// [[Rcpp::export]]
NumericVector MBweights(const IntegerMatrix & ccdata,
			const IntegerVector & ccstatus)
{
  NumericVector rv( ccdata.ncol() );
  unsigned ncontrols = count(ccstatus.begin(),ccstatus.end(),0);
  for( unsigned site = 0 ; site < ccdata.ncol() ; ++site )
    {
      unsigned minor_count = 0;
      for( unsigned ind = 0 ; ind < ccdata.nrow() ; ++ind )
	{
	  if( ccstatus[ind] == 0 )//control
	    {
	      minor_count += ccdata(ind,site);
	    }
	}
      double qi = double(minor_count + 1)/(2.*double(ncontrols)+2.);
      rv[site] = sqrt(double(ccdata.nrow())*qi*(1.-qi) );
    }
  return rv;
}

//' Madsen-Browning test statistics
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return The M-B test statistic for the "general genetic", "recessive", and "dominant" models.
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
//' @details When calculating the rank of an individual's score, the function uses the equivalent of ties="min" in R's rank() function.
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' mbstats = MBstat( rec.ccdata$genos[,which(keep==1)], status )
// [[Rcpp::export]]
Rcpp::List MBstat( const IntegerMatrix & ccdata,
		   const IntegerVector & ccstatus )
{
  stat_MadsenBrowning mb(ccdata.nrow(),count(ccstatus.begin(),ccstatus.end(),0),&ccstatus);
  List rv = stat_calculator(ccdata,ccstatus,mb);
  return rv;
}

//' Get permutation distribution of Madsen-Browning test statistics
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.                                                                                    
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case. 
//' @param nperms The number of permutations to perform
//' @return A data frame of permuted statistics
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' mbstats = MBstat( rec.ccdata$genos[,which(keep==1)], status )
//' mbstats.perm = MB_perm( rec.ccdata$genos[,which(keep==1)], status, 100 )
// [[Rcpp::export]]
DataFrame MB_perm( const IntegerMatrix & ccdata,
		   const IntegerVector & ccstatus,
		   const unsigned & nperms )
{
  NumericVector g(nperms),d(nperms),r(nperms);
  RNGScope scope;
  IntegerVector status = clone(ccstatus);
 
  for( unsigned i = 0 ; i < nperms ; ++i )
    {
      random_shuffle(status.begin(),status.end(),randWrapper);
      List mbstats = MBstat( ccdata, status );
      g[i] = as<double>(mbstats["general"]);
      d[i] = as<double>(mbstats["dominant"]);
      r[i] = as<double>(mbstats["recessive"]);
    }

  return DataFrame::create( Named("general") = g,
			    Named("recessive") = r,
			    Named("dominant") = d );
}
