#ifndef __ESM_HPP__
#define __ESM_HPP__

#include <Rcpp.h>

/*
Association stat from Thornton, Foran, and Long (2013) PLoS Genetics

scores is assumed to be a vector of -log10(p values) for some test of interest

K is the max # of scores to calculate ESM.
 */
double esm( const Rcpp::NumericVector & scores, const unsigned & K );

#endif
