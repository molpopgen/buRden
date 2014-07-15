#ifndef __ESM_FILTER_SITES_HPP__
#define __ESM_FILTER_SITES_HPP__

#include <readCC.hpp>

/*
  Assumes the ccdata is read in with rotate = true (rows are sites, columns are individual
 */
std::vector<short> filter_sites(const CCblock & ccdata,const double & minfreq,const double & maxfreq,const double & rsq_cutoff);

#endif
