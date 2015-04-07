#include <stat_MadsenBrowning.hpp>
#include <functional>
#include <numeric>
#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace std;

stat_MadsenBrowning::stat_MadsenBrowning( const unsigned & __nrows,
					  const unsigned & __ncontrols,
					  const Rcpp::IntegerVector & __ccstatus ) : ncontrols(__ncontrols),
										     minor_count(0),
										     rec_count(0),
										     dom_count(0),
										     ind(0),
										     status(__ccstatus),
										     scores(NumericVector(__nrows)),
										     scores_rec(NumericVector(__nrows)),
										     scores_dom(NumericVector(__nrows)),
										     scores_site(NumericVector(__nrows)),
										     scores_rec_site(NumericVector(__nrows)),
										     scores_dom_site(NumericVector(__nrows))
									   
{
}

void stat_MadsenBrowning::update()
{
  double qi = double(minor_count + 1)/(2.*double(ncontrols)+2.),
    qi_rec = double(rec_count + 1)/(2.*double(ncontrols)+2.),
    qi_dom = double(dom_count + 1)/(2.*double(ncontrols)+2.);
  
  transform( scores_site.begin(),scores_site.begin(),scores_site.begin(),
	     bind2nd( divides<double>(), sqrt(double(scores_site.size())*qi*(1.-qi))) );
  transform( scores_rec_site.begin(),scores_rec_site.begin(),scores_rec_site.begin(),
	     bind2nd( divides<double>(), sqrt(double(scores_rec_site.size())*qi_rec*(1.-qi_rec))) );
  transform( scores_dom_site.begin(),scores_dom_site.begin(),scores_dom_site.begin(),
	     bind2nd( divides<double>(), sqrt(double(scores_dom_site.size())*qi_dom*(1.-qi_dom))) );
  
  scores += scores_site;
  scores_rec += scores_rec_site;
  scores_dom += scores_dom_site;

  //reset variables
  scores_site = NumericVector(scores_site.size());
  scores_rec_site = NumericVector(scores_rec_site.size());
  scores_dom_site = NumericVector(scores_dom_site.size());
  minor_count = rec_count = dom_count = ind = 0;
}

void stat_MadsenBrowning::operator()(const int & genotype,
				     const int & ccstatus)
{
  if( ccstatus == 0 )//control
    {
      switch (genotype)
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
  switch (genotype)
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
  ++ind;
}

Rcpp::List stat_MadsenBrowning::values()
{
  std::vector<double> scores_sorted(scores.begin(),scores.end()),
    scores_rec_sorted(scores_rec.begin(),scores_rec.end()),
    scores_dom_sorted(scores_dom.begin(),scores_dom.end());
  
  sort(scores_sorted.begin(),scores_sorted.end());
  sort(scores_rec_sorted.begin(),scores_rec_sorted.end());
  sort(scores_dom_sorted.begin(),scores_dom_sorted.end());


  double stat=0.,stat_rec=0.,stat_dom=0.;
  for( R_len_t i = 0 ; i < scores.size() ; ++i )
    {
      if( status[i] == 1 )//case
      {
	vector<double>::const_iterator itr = lower_bound(scores_sorted.begin(),scores_sorted.end(),scores[i]);
	//NumericVector::const_iterator itr = find(scores_sorted.begin(),scores_sorted.end(),scores[i]);
	stat += double( itr - scores_sorted.begin() ) + 1.;   //This is the rank
	
	itr = lower_bound(scores_rec_sorted.begin(),scores_rec_sorted.end(),scores_rec[i]);
	stat_rec += double( itr - scores_rec_sorted.begin() ) + 1.;   //This is the rank
	
	itr = lower_bound(scores_dom_sorted.begin(),scores_dom_sorted.end(),scores_dom[i]);
	stat_dom += double( itr - scores_dom_sorted.begin() ) + 1.;   //This is the rank
      }
    }
 
  return List::create( Named("general") = stat,
		       Named("recessive") = stat_rec,
		       Named("dominant") = stat_dom );
}
