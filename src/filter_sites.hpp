#ifndef __FILTER_SITES_HPP__
#define __FILTER_SITES_HPP__

#include <Rcpp.h>

Rcpp::IntegerVector filter_sites(const Rcpp::IntegerMatrix & ccdata,
				 const Rcpp::IntegerVector & ccstatus,
				 const double & minfreq,
				 const double & maxfreq,
				 const double & rsq_cutoff);

#endif
