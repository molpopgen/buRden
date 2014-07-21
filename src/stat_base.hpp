#ifndef __STAT_BASE_HPP__
#define __STAT_BASE_HPP__

#include <Rcpp.h>

class stat_base
{
public:
  virtual void update() = 0;
  virtual void operator()(const int & genotype,
			  const int & ccstatus) = 0;
  virtual Rcpp::List values() = 0;
};

#endif
