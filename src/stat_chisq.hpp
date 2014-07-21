#ifndef __STAT_CHISQ_HPP__
#define __STAT_CHISQ_HPP__

#include <stat_base.hpp>

class stat_chisq : public stat_base
{
  private:
  double log10chisq();
  Rcpp::NumericVector csqs;
  unsigned ctable[4];
public:
  stat_chisq();
  virtual void update();
  virtual void operator()(const int & genotype,
			  const int & ccstatus);
  virtual Rcpp::List values();
};

#endif
