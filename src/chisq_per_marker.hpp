#ifndef __CHISQ_PER_MARKER_HPP__
#define __CHISQ_PER_MARKER_HPP__

#include <readCC.hpp>

/*
  assumes ccdata read in with rotate = true
  so that rows = sites, columns = individuals
 */
std::vector<double> chisq_per_marker( const CCblock * ccdata, 
				      const std::vector<short> * keep, 
				      const std::vector<short>  & ccstatus );

#endif
