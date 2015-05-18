#ifndef __STAT_MADSEN_BROWNING__
#define __STAT_MADSEN_BROWNING__

#include <stat_base.hpp>
#include <vector>
class stat_MadsenBrowning : public stat_base
{
private:
  unsigned ncontrols,minor_count,ind;
  const Rcpp::IntegerVector * status;
  std::vector<double> scores,scores_rec,scores_dom,scores_site,scores_rec_site,scores_dom_site;
public:
  stat_MadsenBrowning( const unsigned & __nrows,
		       const unsigned & __ncontrols,
		       const Rcpp::IntegerVector * __ccstatus); 
  void update();
  void operator()(const int & genotype,
		  const int & ccstatus);
  Rcpp::List values();
};


#endif
