#ifndef __STAT_BASE_HPP__
#define __STAT_BASE_HPP__

#include <Rcpp.h>

/*
  Virtual base class for how to calculate a statistic.
  Tee stat_calculator.cc to make this make a lot more sense.
 */
class stat_base
{
public:
  /*
    After a site is processed, update() is called to update any private
    data that a derived class may contain
   */
  virtual void update() = 0;
  /*
    operator() takes the genotype value for an individual at a certain site,
    and that individuals case/control label as arguments.  It uses
    these arguments to update whatever is necessary to calculate the desired
    statistic
   */
  virtual void operator()(const int & genotype,
			  const int & ccstatus) = 0;
  /*
    values performs any final required calculations
    and returns something interesting
   */
  virtual Rcpp::List values() = 0;
};

#endif
