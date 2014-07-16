#include <ProductMoment.hpp>

using namespace std;
using namespace Rcpp;

//' Pearson's product-moment correlation
//' @param x A vector of values.
//' @param y A vector of values.
//' @return The correlation coefficient between x and y. This should/will be equal to R's cor(x,y).
//' @note This implementation is based on a copy of [libsequence's](http://github.com/molpopgen/libsequence) template function object Sequence::ProductMoment in <Sequence/Correlations.hpp>
// [[Rcpp::export]]
std::iterator_traits<NumericVector::const_iterator>::value_type ProductMoment( const NumericVector & x,
									       const NumericVector & y )
  {
    NumericVector::const_iterator beg_x = x.begin(),
      end_x = x.end(),
      beg_y = y.begin();
    typedef iterator_traits<NumericVector::const_iterator>::value_type rtype;

    if (beg_x >= end_x) return std::numeric_limits<rtype>::min();

    unsigned nsam = 0;
    rtype x_bar=0.;
    rtype y_bar=0.;
    rtype x_squared=0.;
    rtype y_squared=0.;
    rtype xy=0.;
    
    for ( ; beg_x != end_x  ; ++beg_x,++beg_y)
      {
	x_bar += *beg_x;
	x_squared += (*beg_x)*(*beg_x);
	y_bar += *beg_y;
	y_squared += (*beg_y)*(*beg_y);	
	xy += (*beg_x)*(*beg_y);
	++nsam;
      }
    x_bar /= rtype(nsam);
    y_bar /= rtype(nsam);
    x_squared /= rtype(nsam);
    y_squared /= rtype(nsam);
    xy -= x_bar*rtype(nsam)*y_bar;
    xy -= y_bar*rtype(nsam)*x_bar;
    xy += rtype(nsam)*x_bar*y_bar;
    rtype sum_x_sq = rtype(nsam)*(x_squared-x_bar*x_bar);
    rtype sum_y_sq = rtype(nsam)*(y_squared-y_bar*y_bar);
    return(xy/std::pow((sum_x_sq*sum_y_sq),0.5));
  }
