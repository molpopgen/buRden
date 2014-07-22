#include <stat_LiLeal.hpp>
#include <stat_calculator.hpp>
#include <randWrapper.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace std;
stat_LLcollapse::stat_LLcollapse(const double & maf,
				 const Rcpp::IntegerVector & ccstatus,
				 const bool & maf_control) : maf_cutoff(maf),
							     mafc(maf_control),
							     hasRare(IntegerVector(ccstatus.size(),0)),
							     hasRare_site(IntegerVector(ccstatus.size(),0)),
							     status(ccstatus),
							     sum(0),
							     ind(0),
							     ncontrols(count(ccstatus.begin(),ccstatus.end(),0)),
							     N(ccstatus.size())
{
}

void stat_LLcollapse::update()
{
  double maf = double(sum)/double( mafc ? 2*ncontrols : 2*N);

  if ( maf <= maf_cutoff )
    {
      hasRare += hasRare_site;
    }
  sum = ind = 0;
}

void stat_LLcollapse::operator()(const int & genotype,
				 const int & ccstatus)
{
  sum += (mafc) ? (( ccstatus == 0 ) ? genotype : 0) : genotype;
  hasRare_site[ind++] = genotype;
}

Rcpp::List stat_LLcollapse::values()
{
  unsigned co=0,ca=0,cowo=0,cawo=0;
  for( unsigned i = 0 ; i < hasRare.size() ; ++i )
    {
      if( !status[i] ) //control
	{
	  if(hasRare[i])
	    {
	      ++co;
	    }
	  else
	    {
	      ++cowo;
	    }
	}
      else //case
	{
	  if(hasRare[i])
	    {
	      ++ca;
	    }
	  else
	    {
	      ++cawo;
	    }
	}
    }
  double a = co,b=ca,c=cowo,d=cawo;
  double __N = a+b+c+d;
  double rv = log10(__N)+2.*log10(fabs(a*d-b*c)-__N/2.) - ( log10(a+b)+log10(c+d)+log10(b+d)+log10(a+c) );
  return List::create(Named("statistic") = std::pow(10,rv));
}

//' Calculates Li and Leal's collapsed variant statistic, v_c
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param maf Only consider variants whose minor allele frequencies are <= maf
//' @param maf_controls  If true, calculate mafs from controls only.  Otherwise, use all individuals
//' @return A chi-squared statistic based on a 2x2 table of the number of cases and controls with and without rare alleles
//' @references Li, B., & Leal, S. (2008). Methods for detecting associations with rare variants for common diseases: application to analysis of sequence data. The American Journal of Human Genetics, 83(3), 311-321.
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' LL = LLcollapse(rec.ccdata$genos,status,0.01)
// [[Rcpp::export]]
List LLcollapse(const IntegerMatrix & ccdata,
		const IntegerVector & ccstatus,
		const double & maf,
		const bool & maf_controls = false)
{
  stat_LLcollapse f( maf, ccstatus, maf_controls );
  return stat_calculator(ccdata,ccstatus,f);
}

//' Permutation distribution of Li and Leal's collapsed variant statistic, v_c
//' @param ccdata A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param ccstatus A vector of binary phenotype labels.  0 = control, 1 = case.
//' @param nperms The number of permutations to perform
//' @param maf Only consider variants whose minor allele frequencies are <= maf
//' @param maf_controls  If true, calculate mafs from controls only.  Otherwise, use all individuals
//' @return The non-centrality parameter of a chi-squared distribution.  This is obtained using the proportion of controls and cases with rare variants.
//' @references Li, B., & Leal, S. (2008). Methods for detecting associations with rare variants for common diseases: application to analysis of sequence data. The American Journal of Human Genetics, 83(3), 311-321.
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' LL.perm = LLcollapse_perm(rec.ccdata$genos,status,10,0.01)
// [[Rcpp::export]]		  
NumericVector LLcollapse_perm(const IntegerMatrix & ccdata,
			      const IntegerVector & ccstatus,
			      const unsigned & nperms,
			      const double & maf,
			      const bool & maf_controls = false)
{
  NumericVector rv(nperms);
  RNGScope scope;
  IntegerVector status = clone(ccstatus);

  for( unsigned i = 0 ; i < nperms ; ++i )
    {
      random_shuffle(status.begin(),status.end(),randWrapper);
      stat_LLcollapse f( maf, status, maf_controls );
      rv[i] = as<double>( stat_calculator(ccdata,status,f)["statistic"] );
    }
  return rv;
}
