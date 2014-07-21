#ifndef __STAT_MADSEN_BROWNING__
#define __STAT_MADSEN_BROWNING__

#include <stat_base.hpp>

class stat_MadsenBrowning : public stat_base
{
private:
  unsigned ncontrols,minor_count,rec_count,dom_count,ind;
  Rcpp::IntegerVector status;
  Rcpp::NumericVector scores,scores_rec,scores_dom,scores_site,scores_rec_site,scores_dom_site;
public:
  stat_MadsenBrowning( const unsigned & __nrows,
		       const unsigned & __ncontrols,
		       const Rcpp::IntegerVector & __ccstatus); 
  virtual void update();
  virtual void operator()(const int & genotype,
			  const int & ccstatus);
  virtual Rcpp::List values();
};


#endif
