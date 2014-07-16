#ifndef __CHISQ_PER_MARKER_HPP__
#define __CHISQ_PER_MARKER_HPP__

#include <Rcpp.h>

Rcpp::NumericVector chisq_per_marker( const Rcpp::IntegerMatrix & ccdata,
				      const Rcpp::IntegerVector & ccstatus );


#endif
