#ifndef __MBstat_HPP__
#define __MBstat_HPP__

#include <Rcpp.h>

Rcpp::List MBstat( const Rcpp::IntegerMatrix & data,
		   const Rcpp::IntegerVector & status );

#endif
