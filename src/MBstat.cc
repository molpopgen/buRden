#include<MBstat.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

using namespace Rcpp;
using namespace std;

//' Calculate Madsen-Browning weights.
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return An array of weights, one for each column in data.
//' @details Calculation is done under the "general genetic model" defined in Madsen and Browning.
// [[Rcpp::export]]
NumericVector MBweights(const IntegerMatrix & data,
			const IntegerVector & status)
{
  NumericVector rv( data.ncol() );
  unsigned ncontrols = count(status.begin(),status.end(),0);
  for( unsigned site = 0 ; site < data.ncol() ; ++site )
    {
      unsigned minor_count = 0;
      for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
	{
	  if( status[ind] == 0 )//control
	    {
	      minor_count += data(ind,site);
	    }
	}
      double qi = double(minor_count + 1)/(2.*double(ncontrols)+2.);
      rv[site] = sqrt(double(ncontrols)*qi*(1.-qi) );
    }
  return rv;
}

//' Madsen-Browning test statistics
//' @param data A matrix of markers (columns) and individuals (rows).  Data are coded as the number of copies of the minor allele.
//' @param status A vector of binary phenotype labels.  0 = control, 1 = case.
//' @return The M-B test statistic for the "general genetic", "recessive", and "dominant" models.
//' @examples
//' data(rec.ccdata)
//' status = c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases))
//' #filter out common alleles and marker pairs in high LD
//' keep = filter_sites(rec.ccdata$genos,status,0,0.05,0.8)
//' mbstats = MBstat( rec.ccdata$genos[,which(keep==1)], status )
// [[Rcpp::export]]
Rcpp::List MBstat( const IntegerMatrix & data,
		   const IntegerVector & status )
{
  NumericVector weights = MBweights(data,status);
  unsigned ncontrols = count(status.begin(),status.end(),0);
  //Here, we'll keep track of overall scores
  NumericVector scores(data.nrow()),
    scores_rec(data.nrow()),scores_dom(data.nrow());
  for( unsigned site = 0 ; site < data.ncol() ; ++site )
    {
      unsigned minor_count = 0,rec_count=0,dom_count=0;
      NumericVector scores_site(data.nrow()),
	scores_rec_site(data.nrow()),scores_dom_site(data.nrow());
      for( unsigned ind = 0 ; ind < data.nrow() ; ++ind )
	{
	  if( status[ind] == 0 )//control
	    {
	      switch (data(ind,site))
		{
		case 0:
		  break;
		case 1:
		  scores_site[ind]+=1;
		  scores_dom_site[ind]+=1;
		  minor_count+=1;
		  dom_count += 1;
		  break;
		case 2:
		  scores_site[ind]+=2;
		  scores_dom_site[ind]+=1;
		  scores_rec_site[ind]+=1;
		  minor_count += 2;
		  dom_count += 1;
		  rec_count += 1;
		  break;
		default:
		  stop("MBstat: genotype code other than 0, 1, or 2 is not allowed!");
		}
	    }
	  switch (data(ind,site))
	    {
	    case 0:
	      break;
	    case 1:
	      scores_site[ind]+=1;
	      scores_dom_site[ind]+=1;
	      break;
	    case 2:
	      scores_site[ind]+=2;
	      scores_dom_site[ind]+=1;
	      scores_rec_site[ind]+=1;
	      break;
	    default:
	      stop("MBstat: genotype code other than 0, 1, or 2 is not allowed!");
	    }
	}
      double qi = double(minor_count + 1)/(2.*double(ncontrols)+2.),
	qi_rec = double(rec_count + 1)/(2.*double(ncontrols)+2.),
	qi_dom = double(dom_count + 1)/(2.*double(ncontrols)+2.);
  
      transform( scores_site.begin(),scores_site.begin(),scores_site.begin(),
		 bind2nd( divides<double>(), qi ) );
      transform( scores_rec_site.begin(),scores_rec_site.begin(),scores_rec_site.begin(),
		 bind2nd( divides<double>(), qi_rec ) );
      transform( scores_dom_site.begin(),scores_dom_site.begin(),scores_dom_site.begin(),
		 bind2nd( divides<double>(), qi_dom ) );
  
      scores += scores_site;
      scores_rec += scores_rec_site;
      scores_dom += scores_dom_site;
    }

  vector<double> scores_sorted(scores.begin(),scores.end()),
    scores_rec_sorted(scores_rec.begin(),scores_rec.end()),
    scores_dom_sorted(scores_dom.begin(),scores_dom.end());

  sort(scores_sorted.begin(),scores_sorted.end());
  sort(scores_rec_sorted.begin(),scores_rec_sorted.end());
  sort(scores_dom_sorted.begin(),scores_dom_sorted.end());
  
  double stat=0.,stat_rec=0.,stat_dom=0.;
  for( unsigned i = 0 ; i < scores.size() ; ++i )
    {
      if( status[i] == 1 )//control
	{
	  vector<double>::const_iterator itr = find(scores_sorted.begin(),scores_sorted.end(),scores[i]);
	  stat += double( itr - scores_sorted.begin() ) + 1.;   //This is the rank

	  itr = find(scores_rec_sorted.begin(),scores_rec_sorted.end(),scores_rec[i]);
	  stat_rec += double( itr - scores_rec_sorted.begin() ) + 1.;   //This is the rank

	  itr = find(scores_dom_sorted.begin(),scores_dom_sorted.end(),scores_dom[i]);
	  stat_dom += double( itr - scores_dom_sorted.begin() ) + 1.;   //This is the rank
	}
    }
  return List::create( Named("general") = stat,
		       Named("recessive") = stat_rec,
		       Named("dominant") = stat_dom );
}
