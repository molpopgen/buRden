#include <randWrapper.hpp>

#include <Rcpp.h>
#include <cmath>
int randWrapper(const int & n) { return floor(R::runif(0.,1.)*n); }
