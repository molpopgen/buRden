#ifndef __STAT_CALPHA_HPP__
#define __STAT_CALPHA_HPP__

#include <stat_base.hpp>
#include <map>
class stat_cAlpha : public stat_base
{
private:
  double T,p0;
  unsigned n_i,y_i;
  bool norm,simple;
  std::map<unsigned,unsigned> ns;
  double Z() const;
public:
  stat_cAlpha(const Rcpp::IntegerVector & status,
	      const bool & normalize = false,
	      const bool & simplecounts = false);
  virtual void update();
  virtual void operator()(const int & genotype,
			  const int & ccstatus);
  virtual Rcpp::List values();
};

#endif
