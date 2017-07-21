#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List cppswitches(int k,int NUMA,int maxdonors,bool THIN,int NUMP,int G,IntegerVector NL,IntegerVector label,NumericVector sumfors,
    NumericVector backs,NumericVector transitions,IntegerVector flips,NumericVector mutmat,int maxmiss,List UMATCH,int max_umatch_size,List dw, List tw,
    IntegerVector gobs,IntegerVector kndonors,IntegerVector donates)
{
  int offset=(THIN ? maxdonors : 0); // setting offset=0 for maxdonors==NUMP => use single vector of donors
  maxmiss=maxmiss+1; // the +1 is to make room for the zero
  IntegerVector dim(2);dim[0]=G;dim[1]=NUMP;
  Dimension d(dim);
  NumericVector unscaledswitches(d); // unscaled probs definitely switching i.e. to other hap
  NumericVector switches(d); // scaled probs definitely switching i.e. to other hap
  IntegerVector dims(1);dims[0]=NUMP;
  Dimension ds(dims);
  double invsum, tmp;
  int g,jl,jk,djk,gum,glk,dg,ndonors=maxdonors;
  int h=k, otherh=((k%2) ? k+1 : k-1);
  k=k-1;otherh=otherh-1; // move to c++ indexing
  List ldw_w=as<List>(dw[0]), ltw_w=as<List>(tw[0]);
  IntegerVector dw_w(NUMP); // many g will be blank, if obs then length is NUMP
  IntegerVector tw_w(NUMA); // many g will be blank, if obs then length is NUMA
  IntegerVector umatch(max_umatch_size); // many g will be blank, if obs then length is NUMA
  IntegerVector dw_size=as<IntegerVector>(dw[1]);
  for (jk=0;jk<NUMP;jk++)
    label[jk]=label[jk]-1; // re-index to cpp
  for (g=1; g<G; g++) {
    dw_w=as<IntegerVector>(ldw_w[g]);
    tw_w=as<IntegerVector>(ltw_w[g]);
    umatch=as<IntegerVector>(UMATCH[g]);
    h=(flips[g] ? otherh : k); // if flips[g]==0 then h=same as was else with even k to k-1 and odd k to k+1
    ndonors=kndonors[g];
    invsum=0.0;
    for (djk=0;djk<ndonors;djk++) {
      dg=g*offset+djk;
      jk=donates[dg]; 
      jl=label[jk];
      tmp=sumfors[(g-1)]*backs[g*maxdonors+djk]*transitions[jl*2];
      if (gobs[g]>0) 
      {
	gum=(tw_w[h])*dw_size[g]+dw_w[jk];
	glk=umatch[gum];
        tmp*=mutmat[maxmiss*glk+(gobs[g]-glk)]; 
      }
      unscaledswitches[g*NUMP+jk]=tmp;
      invsum+=tmp;
    }
    invsum=1.0/invsum;
    for (jk=0;jk<NUMP;jk++)
      switches[g*NUMP+jk]+=unscaledswitches[g*NUMP+jk]*invsum;
  }
  List ret;
  ret["switches"]=switches;
  return wrap(ret);
}
