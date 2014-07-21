#include<MBstat.hpp>
#include <stat_MadsenBrowning.hpp>
#include <stat_calculator.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

using namespace Rcpp;
using namespace std;

//' Calculate Madsen-Browning weights.
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return An array of weights, one for each column in data.
//' @details Calculation is done under the "general genetic model" defined in Madsen and Browning.
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
// [[Rcpp::export]]
NumericVector MBweights(const IntegerMatrix & data,
			const IntegerVector & status)
{
  NumericVector rv( data.ncol() );
  unsigned ncontrols = count(status.begin(),status.end(),0);
  for( unsigned site = 0 ; site < data.ncol() ; ++site )
    {
      unsigned minor_count = 0;
      for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
	{
	  if( status[ind] == 0 )//control
	    {
	      minor_count += data(ind,site);
	    }
	}
      double qi = double(minor_count + 1)/(2.*double(ncontrols)+2.);
      rv[site] = sqrt(double(data.nrow())*qi*(1.-qi) );
    }
  return rv;
}

//' Madsen-Browning test statistics
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
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
Rcpp::List MBstat( const IntegerMatrix & data,
		   const IntegerVector & status )
{
  stat_MadsenBrowning mb(data.nrow(),count(status.begin(),status.end(),0),status);
  List rv = stat_calculator(data,status,mb);
  return rv;
}
