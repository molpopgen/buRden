#ifndef __STAT_CALCULATOR_HPP__
#define __STAT_CALCULATOR_HPP__

#include <Rcpp.h>
#include <stat_base.hpp>

Rcpp::List stat_calculator(const Rcpp::IntegerMatrix & data,
			   const Rcpp::IntegerVector & status,
			   stat_base & f);

#endif
