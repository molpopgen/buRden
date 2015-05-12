#include <esm.hpp>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>

#include <boost/tuple/tuple.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/tail.hpp>
#include <boost/accumulators/statistics/tail_variate.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
using namespace std;
//using namespace Rcpp;
using namespace boost::accumulators;

namespace {
  /*
    Awaiting the days when CRAN will take cpp11...
  */
  struct __esm_accumulator {
    mutable unsigned i,ntests;
    __esm_accumulator(const unsigned & n) : i(0u),ntests(n) {}
    typedef double result_type;
    inline result_type operator()(const double & a, const double & b) const
    {
      double rv = a + (b - (-log10((i+1)/double(ntests))));
      ++i;
      return rv;
    }
  };
}

namespace {
  /*
    Awaiting the days when CRAN will take cpp11...
  */
  struct __esmw_accumulator {
    mutable unsigned i,ntests;
    __esmw_accumulator(const unsigned & n) : i(0u),ntests(n) {}
    typedef double result_type;
    inline result_type operator()(const double & a, const boost::tuple<const double &, const double & > & b) const
    {
      
      double rv = a + (get<0>(b) + log10((i+1)/double(ntests)))*get<1>(b);
      ++i;
      return rv;
    }
  };
}

//' Association stat from Thornton, Foran, and Long (2013) PLoS Genetics
//' @param scores A vector of single-marker association test scores, on a -log10 scale
//' @param K the number of markers used to calculate ESM_K
//' @return The ESM_K test statistic value
//' @references Thornton, K. R., Foran, A. J., & Long, A. D. (2013). Properties and Modeling of GWAS when Complex Disease Risk Is Due to Non-Complementing, Deleterious Mutations in Genes of Large Effect. PLoS Genetics, 9(2), e1003258. doi:10.1371/journal.pgen.1003258
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' rec.ccdata.chisq = chisq_per_marker(rec.ccdata$genos[,which(keep==1)],status)
//' rec.ccdata.esm = esm( rec.ccdata.chisq, 50 )
// [[Rcpp::export]]
double esm( const Rcpp::NumericVector & scores, const unsigned & K )
{
  typedef tag::tail_variate< int, tag::covariate1, boost::accumulators::right > my_tail_variate_tag;
  unsigned ntests = scores.size();
  accumulator_set<double, stats<tag::tail<boost::accumulators::right> > > acc(tag::tail<boost::accumulators::right>::cache_size =std::min(K,ntests));
  for( Rcpp::NumericVector::const_iterator itr = scores.begin() ; 
       itr != scores.end() ; ++itr )
    {
      acc(*itr);
    }
  __esm_accumulator b(ntests);
  double rv = accumulate(tail(acc).begin(),tail(acc).end(),0.,b);
  return rv;
}

//' weighted verstion of Association stat from Thornton, Foran, and Long (2013) PLoS Genetics
//' @param scores A vector of single-marker association test scores, on a -log10 scale
//' @param weights A vector of weights to use for each marker, such as dbetat(MAF,1,25)
//' @param K the number of markers used to calculate ESM_K
//' @return The ESM_K test statistic value
//' @references Thornton, K. R., Foran, A. J., & Long, A. D. (2013). Properties and Modeling of GWAS when Complex Disease Risk Is Due to Non-Complementing, Deleterious Mutations in Genes of Large Effect. PLoS Genetics, 9(2), e1003258. doi:10.1371/journal.pgen.1003258
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' MAF = colSums(ccdata$genos[,which(keep==1)])/(2*nrow(ccdata$genos))
//' b.weight = dbeta(MAF,1,25)
//' lin.weight = rep(1,length(MAF))
//' rec.ccdata.chisq = chisq_per_marker(rec.ccdata$genos[,which(keep==1)],status)
//' rec.ccdata.esm = esm( rec.ccdata.chisq,50 )
//' rec.ccdata.esm.lin =  esmw( rec.ccdata.chisq,lin.weight, 50 )
//' rec.ccdata.esm.beta = esmw( rec.ccdata.chisq,b.weight, 50 )
// [[Rcpp::export]]
double esmw( const Rcpp::NumericVector & scores, const Rcpp::NumericVector & weights, const unsigned & K )
{
  typedef tag::tail_variate< boost::tuple<double,double>, tag::covariate1, boost::accumulators::right > my_tail_variate_tag;
  unsigned ntests = scores.size();
  accumulator_set<double, features< my_tail_variate_tag > > acc(tag::tail<boost::accumulators::right>::cache_size =std::min(K,ntests));
  Rcpp::NumericVector::const_iterator itr1;
  Rcpp::NumericVector::const_iterator itr2;
  for(itr1 = scores.begin(),itr2 = weights.begin() ; (itr1 != scores.end()) && (itr2 !=weights.end()) ; ++itr1, ++itr2  )
    {
      acc(*itr1,covariate1= boost::make_tuple(*itr1,*itr2));
    }
  extractor< my_tail_variate_tag > my_tail_variate;
  __esmw_accumulator b(ntests);
  double rv = accumulate(my_tail_variate(acc).begin(),my_tail_variate(acc).end(),0.,b);
  return rv;
}
