#include <Rcpp.h>
#include <stat_allstats.hpp>
#include <stat_calculator.hpp>
#include <randWrapper.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace std;

//' Calculate all burden statistics simultaneously
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param esm_K The number of markers to use in the calculation of ESM_K
//' @param normalize_calpha If TRUE, return T/sqrt(Z), otherwise return T.
//' @param simplecount_calpha see Details
//' @return A list of values for all burden statistics
//' @references Neale, B. M., Rivas, M. A., Voight, B. F., Altshuler, D., Devlin, B., Orho-Melander, M., et al. (2011). Testing for an Unusual Distribution of Rare Variants. PLoS Genetics, 7(3), e1001322. doi:10.1371/journal.pgen.1001322
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
//' @references Thornton, K. R., Foran, A. J., & Long, A. D. (2013). Properties and Modeling of GWAS when Complex Disease Risk Is Due to Non-Complementing, Deleterious Mutations in Genes of Large Effect. PLoS Genetics, 9(2), e1003258. doi:10.1371/journal.pgen.1003258
//' @details  When simplecount_alpha = FALSE, heterozygous and homozygous genotypes are treated as different numbers of observations
//' of the mutation.  In other wordes, simplecounts = FALSE is equivalent to colSums( ccdata[status==1,] ).  When simplecounts=TRUE,
//' all nonzero genotype values are treated as the value 1, equivalent to  apply(data[status==1,], 2, function(x) sum(x>0, na.rm=TRUE)).
//' The latter method is used by the R package AssotesteR.
// [[Rcpp::export]]
List allBurdenStats( const IntegerMatrix & ccdata,
		     const IntegerVector & ccstatus,
		     const unsigned & esm_K,
		     const bool normalize_calpha = false,
		     const bool simplecount_calpha = false )
{
  stat_allstats f(ccdata.nrow(),ccstatus,esm_K,normalize_calpha,simplecount_calpha);
  return stat_calculator(ccdata,ccstatus,f);
}

//' Estimate p-values for all burden statistics by permutation
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param esm_K The number of markers to use in the calculation of ESM_K
//' @param nperms Number of permutations to perform
//' @param normalize_calpha If TRUE, return T/sqrt(Z), otherwise return T.
//' @param simplecount_calpha see Details
//' @return A list of p-values for all burden statistics.
//' @references Neale, B. M., Rivas, M. A., Voight, B. F., Altshuler, D., Devlin, B., Orho-Melander, M., et al. (2011). Testing for an Unusual Distribution of Rare Variants. PLoS Genetics, 7(3), e1001322. doi:10.1371/journal.pgen.1001322
//' @references Madsen, B. E., & Browning, S. R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genetics, 5(2), e1000384. doi:10.1371/journal.pgen.1000384
//' @references Thornton, K. R., Foran, A. J., & Long, A. D. (2013). Properties and Modeling of GWAS when Complex Disease Risk Is Due to Non-Complementing, Deleterious Mutations in Genes of Large Effect. PLoS Genetics, 9(2), e1003258. doi:10.1371/journal.pgen.1003258
//' @details  When simplecount_alpha = FALSE, heterozygous and homozygous genotypes are treated as different numbers of observations
//' of the mutation.  In other wordes, simplecounts = FALSE is equivalent to colSums( ccdata[status==1,] ).  When simplecounts=TRUE,
//' all nonzero genotype values are treated as the value 1, equivalent to  apply(data[status==1,], 2, function(x) sum(x>0, na.rm=TRUE)).
//' The latter method is used by the R package AssotesteR.
// [[Rcpp::export]]
List allBurdenStatsPerm( const IntegerMatrix & ccdata,
			 const IntegerVector & ccstatus,
			 const unsigned & esm_K,
			 const unsigned & nperms,
			 const bool normalize_calpha = false,
			 const bool simplecount_calpha = false )
{
  RNGScope scope;
  IntegerVector status = clone(ccstatus);
  //store permutation distributions
  NumericVector esm_p(nperms),
    calpha_p(nperms),
    MBg_p(nperms),
    MBr_p(nperms),
    MBd_p(nperms);

  for( unsigned i = 0 ; i < nperms; ++i )
    {
      random_shuffle(status.begin(),status.end(),randWrapper); 
      stat_allstats f(ccdata.nrow(),status,esm_K,normalize_calpha,simplecount_calpha);
      List perm_vals = stat_calculator(ccdata,status,f);
      esm_p[i] = as<double>( perm_vals["esm.stat"] );
      calpha_p[i] = as<double>( perm_vals["calpha.stat"] );
      MBg_p[i] = as<double>( perm_vals["MB.general.stat"] );
      MBr_p[i] = as<double>( perm_vals["MB.recessive.stat"] );
      MBd_p[i] = as<double>( perm_vals["MB.dominant.stat"] );
    }
  return List::create( Named("esm.permdist") = esm_p,
		       Named("calpha.permdist") = calpha_p,
		       Named("MB.general.permdist") = MBg_p,
		       Named("MB.recessive.permdist") = MBr_p,
		       Named("MB.dominant.permdist") = MBd_p 
		       );
}
