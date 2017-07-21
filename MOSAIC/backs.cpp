#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
void cppbackward(int k,int NUMA,int maxdonors,bool THIN,int NUMP,int L, int gl,int gu,
    int G, NumericVector transitions,List UMATCH,int max_umatch_size,List dw,List tw,IntegerVector gobs,NumericVector mutmat,int maxmiss,
    IntegerVector label,IntegerVector kndonors,IntegerVector donates,IntegerVector donatesr,
    IntegerVector flips,NumericVector backs,NumericVector scalefactor)
{
  int offset=(THIN ? maxdonors : 0); // setting offset=0 for maxdonors==NUMP => use single vector of donors
  maxmiss=maxmiss+1; // the +1 is to make room for the zero
  // end of arguments passed 
  double tmpsum=0.0, tmpsum2=0.0;
  int NNL=NUMP*L;
  int i,l,g,dlk,lk2,lk,llk,ndonors=maxdonors;
  int dg, glk, gum, gm, im,lm, ilk, lkl, il; // convenience indices for faster indexing
  int h=k, otherh=((k%2) ? k+1 : k-1);
  k=k-1;otherh=otherh-1; // move to c++ indexing
  List ldw_w=as<List>(dw[0]), ltw_w=as<List>(tw[0]);
  IntegerVector dw_w(NUMP); // many g will be blank, if obs then length is NUMP
  IntegerVector tw_w(NUMA); // many g will be blank, if obs then length is NUMA
  IntegerVector umatch(max_umatch_size); // many g will be blank, if obs then length is NUMA
  IntegerVector dw_size=as<IntegerVector>(dw[1]);
  IntegerVector dim(2); dim[0]=L; dim[1]=L;
  Dimension d(dim);
  NumericVector sumstates(d);
  NumericVector emissions(NNL);
  g=gu-1;
  for (lk=0;lk<NUMP;lk++)
    label[lk]=label[lk]-1;
  if (gu==G) {
    gm=g*maxdonors*L;
    ndonors=kndonors[g];
    for (l=0;l<L;l++) {
      lm=l*maxdonors;
      for (dlk=0;dlk<ndonors;dlk++) {
	backs[gm+lm+dlk]=1.0;
      }
    }
    scalefactor[g]=1.0/NNL;
  }

  for (g=gu-1;g>gl;g--) {
    dw_w=as<IntegerVector>(ldw_w[g]); 
    tw_w=as<IntegerVector>(ltw_w[g]); 
    umatch=as<IntegerVector>(UMATCH[g]); 
    gm=g*maxdonors*L;
    ndonors=kndonors[g];
    h=(flips[g] ? otherh : k); // if flips[g]==0 then h=same as was else with even k to k-1 and odd k to k+1
    for (i=0;i<L;i++) {
      lm=i*maxdonors;
      im=i*NUMP;
      for (l=0;l<L;l++)
	sumstates[l*L+i]=0.0;
      for (dlk=0;dlk<ndonors;dlk++) {
	dg=g*offset+dlk; 
	lk=donates[dg];
	llk=label[lk];
  	il=im+lk;
  	lkl=2*L*(llk*L + i);
	if (gobs[g]==0) 
	{
	  emissions[il]=1.0;
	  for (l=0;l<L;l++) {
	    ilk=lkl + l; // transition index
	    sumstates[l*L+i]+=backs[gm+lm+dlk]*transitions[ilk];
	  }
	} else 
	{
	  gum=(tw_w[h])*dw_size[g]+dw_w[lk]; 
	  glk=umatch[gum]; 
	  emissions[il]=mutmat[maxmiss*L*glk+L*(gobs[g]-glk)+i]; // i appears here in case we want anc-specific theta
	  for (l=0;l<L;l++) {
	    ilk=lkl + l; // transition index
	    sumstates[l*L+i]+=backs[gm+lm+dlk]*transitions[ilk]*emissions[il];
	  }
	}
      }
    }
    ndonors=kndonors[g-1];
    tmpsum=0.0;
    for (dlk=0;dlk<ndonors;dlk++) {
      dg=(g-1)*offset+dlk; 
      lk=donates[dg];
      llk=label[lk];
      lk2=donatesr[dg];
      for (l=0;l<L;l++) {
	tmpsum2=0.0; // one per l and lk
	for (i=0;i<L;i++) {
	  lm=i*maxdonors;
	  im=i*NUMP;
	  il=im+lk;
	  tmpsum2+=sumstates[l*L+i]; // add in all switch terms
	  if (lk2!=-1) 
	  {
	    lkl=2*L*(llk*L + i);
	    if (emissions[il]==0) emissions[il]=1.0; // this can happen if more donors here but no obs
	    ilk=lkl + l; // transition index
	    tmpsum2+=backs[gm+lm+lk2]*(transitions[ilk+L]-transitions[ilk])*emissions[il]; // remove wrong switch term
	  }
	}
	//if (tmpsum2<1.0e-16) tmpsum2=1.0e-16;
	backs[gm-maxdonors*L+l*maxdonors+dlk]=tmpsum2; // -maxdonors*L takes us to g-1
	tmpsum+=tmpsum2;
      }
    }
    scalefactor[g-1]=1.0/tmpsum;
    for (l=0;l<L;l++) {
      lm=l*maxdonors;
      for (dlk=0;dlk<ndonors;dlk++) {
	backs[gm-maxdonors*L+lm+dlk]*=scalefactor[g-1]; // -maxdonors*L takes us to g-1
      }
    }
  }
}
