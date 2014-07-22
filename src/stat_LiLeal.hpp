#ifndef __STAT_LI_LEAL_HPP__
#define __STAT_LI_LEAL_HPP__

#include <stat_base.hpp>

/*
  Their notation:
  N = Number of individuals
  \phi_A = Proportion of cases carrying a rare variant
  \phi_A^- = 
  
 */

//Collapsing method on p. 314
class stat_LLcollapse : public stat_base
{
private:
  mutable double maf_cutoff;
  mutable bool mafc;
  Rcpp::IntegerVector hasRare,hasRare_site,status;
  unsigned sum,ind,ncontrols,N;
public:
  stat_LLcollapse( const double & __maf, const Rcpp::IntegerVector & ccstatus,
		   const bool & maf_control = true );
  virtual void update();
  virtual void operator()(const int & genotype,
			  const int & ccstatus);
  virtual Rcpp::List values();
};

#endif
