#ifndef __ProductMoment_HPP__
#define __ProductMoment_HPP__
#include <Rcpp.h>
#include <iterator>

std::iterator_traits<Rcpp::NumericVector::const_iterator>::value_type ProductMoment( const Rcpp::NumericVector & x,
										     const Rcpp::NumericVector & y );

#endif
