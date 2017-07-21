#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List cppprobs(int k,int NUMA,int maxdonors,bool THIN,int L,int kLL,int NN,int NUMP,int G,IntegerVector label,NumericVector fors,
    NumericVector sumfors,NumericVector backs,NumericVector transitions,IntegerVector flips,NumericVector mutmat,int maxmiss,
    List UMATCH, int max_umatch_size, List dw, List tw, IntegerVector gobs,IntegerVector kndonors,IntegerVector donates,IntegerVector donatesl)
{
  int offset=(THIN ? maxdonors : 0); // setting offset=0 for maxdonors==NUMP => use single vector of donors
  maxmiss=maxmiss+1; // the +1 is to make room for the zero
  IntegerVector dim(3);dim[0]=L;dim[1]=kLL;dim[2]=L;
  Dimension d(dim);
  NumericVector unscaledswitches(d); // unscaled probs definitely switching i.e. to other hap
  NumericVector switches(d); // scaled probs definitely switching i.e. to other hap
  NumericVector errors(L), tmperrors(L); // counts of errors in copying
  NumericVector hits(L), tmphits(L); // counts of hits in copying
  IntegerVector dims(2);dims[0]=NUMP;dims[1]=L;
  Dimension ds(dims);
  NumericVector self(ds); // scaled probs perhaps switching i.e. same hap and anc
  NumericVector unscaledself(ds); // unscaled probs perhaps switching i.e. same hap and anc
  double invsum, tmp, ginvsum;
  int g,ia,ja,jl,jk,djk,jk2,gum,glk,dg,ndonors=maxdonors;
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
    ginvsum=0.0;
    invsum=0.0;
    for (ia=0;ia<L;ia++) {
      tmperrors[ia]=0.0;tmphits[ia]=0.0;
      for (ja=0;ja<L;ja++) {
	for (jl=0;jl<kLL;jl++)
	  unscaledswitches[ja*kLL*L + jl*L + ia]=0.0;
	for (djk=0;djk<ndonors;djk++) {
	  dg=g*offset+djk;
	  jk=donates[dg]; 
	  jl=label[jk];
	  tmp=sumfors[(g-1)*L+ia]*backs[g*maxdonors*L+ja*maxdonors+djk]*
	    transitions[jl*L*2*L + ja*2*L + 0*L + ia]; // switches only; remove switch-to-self from this below
	  if (gobs[g]>0) 
	  {
	    gum=(tw_w[h])*dw_size[g]+dw_w[jk];
	    glk=umatch[gum]; 
	    tmp*=mutmat[maxmiss*L*glk+L*(gobs[g]-glk)+ja];
	  }
	  unscaledswitches[ja*kLL*L + jl*L + ia]+=tmp;
	  invsum+=tmp;
	}
      }
      for (djk=0;djk<ndonors;djk++) {
	dg=g*offset+djk;
	jk=donates[dg];
	jl=label[jk];
	jk2=donatesl[dg];
	tmp=0.0;
	if (jk2!=-1) 
	  tmp=fors[(g-1)*maxdonors*L+ia*maxdonors+jk2];
	tmp*=backs[g*maxdonors*L+ia*maxdonors+djk];
	if (gobs[g]>0) 
	{
	  gum=(tw_w[h])*dw_size[g]+dw_w[jk];
	  glk=umatch[gum]; 
	  tmp*=mutmat[maxmiss*L*glk+L*(gobs[g]-glk)+ia];
	}
	unscaledswitches[ia*kLL*L + jl*L + ia]-=tmp*transitions[jl*L*2*L + ia*2*L + 0*L + ia]; // remove switch-to-self from switches above to avoid double counting
	invsum-=tmp*transitions[jl*L*2*L + ia*2*L + 0*L + ia]; // and remove switch-to-self from running total count
	unscaledself[ia*NUMP+jk]=tmp*transitions[jl*L*2*L + ia*2*L + 1*L + ia]; // includes both switches-to-self and non-switches
	invsum+=unscaledself[ia*NUMP + jk]; 
	tmp=fors[g*maxdonors*L+ia*maxdonors+djk]*backs[g*maxdonors*L+ia*maxdonors+djk];
	ginvsum+=tmp;
	if (gobs[g]>0)
	{
  	  tmperrors[ia]+=tmp*(gobs[g]-glk);
	  tmphits[ia]+=tmp*glk;
	}
      }
    }
    invsum=1.0/invsum;
    ginvsum=1.0/ginvsum;
    for (ia=0;ia<L;ia++)
    {
      for (ja=0;ja<L;ja++)
	for (jl=0;jl<kLL;jl++)
	  switches[ja*kLL*L + jl*L + ia]+=unscaledswitches[ja*kLL*L + jl*L + ia]*invsum;
      for (djk=0;djk<ndonors;djk++) {
	dg=g*offset+djk;
	jk=donates[dg]; 
	self[ia*NUMP + jk]+=unscaledself[ia*NUMP + jk]*invsum;
      }
      tmperrors[ia]*=ginvsum;
      tmphits[ia]*=ginvsum;
      errors[ia]+=tmperrors[ia];
      hits[ia]+=tmphits[ia];
    }
  }
  List ret;
  ret["switches"]=switches;
  ret["self"]=self;
  ret["errors"]=errors;
  ret["hits"]=hits;
  return wrap(ret);
}
