#include <Rcpp.h>
#include <esm.hpp>
#include <chisq_per_marker.hpp>

using namespace Rcpp;

//' Association stat from Thornton, Foran, and Long (2013) PLoS Genetics
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param k The number of markers for the ESM_K statistic.
//' @return The ESM_K test statistic value based on chi-squared tests per marker.
//' @references Thornton, K. R., Foran, A. J., & Long, A. D. (2013). Properties and Modeling of GWAS when Complex Disease Risk Is Due to Non-Complementing, Deleterious Mutations in Genes of Large Effect. PLoS Genetics, 9(2), e1003258. doi:10.1371/journal.pgen.1003258
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' rec.ccdata.chisq = esm_chisq( rec.ccdata$genos[,which(keep==1)], status, 50 )
// [[Rcpp::export]]
double esm_chisq( const IntegerMatrix & ccdata,
		  const IntegerVector & ccstatus,
		  const unsigned & k)
{
  NumericVector c = chisq_per_marker( ccdata, ccstatus );
  double stat = esm(c,k);
  return( stat );
}
