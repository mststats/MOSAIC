#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
void cppforward(int k,int NUMA,int maxdonors,bool THIN,int NUMP,int kLL,int L, int gl,int gu,
    int G, NumericVector transitions,List UMATCH,int max_umatch_size,List dw,List tw,IntegerVector gobs,NumericVector mutmat,int maxmiss,
    NumericVector initProb,IntegerVector label,IntegerVector kndonors,IntegerVector donates,IntegerVector donatesl,
    IntegerVector flips,NumericVector fors,NumericVector sumfors,NumericVector scalefactor)
{
  int offset=(THIN ? maxdonors : 0); // setting offset=0 for maxdonors==NUMP => use single vector of donors
  maxmiss=maxmiss+1; // the +1 is to make room for the zero
  // end of arguments passed 
  double tmpsum=0.0, tmpsum2=0.0;
  int i,l,g,dlk,lk2,lk,llk,ndonors=maxdonors; 
  int dg, gum, glk, gm, lm, gld, ilk, lkl; // convenience indices for faster indexing
  int h=k, otherh=((k%2) ? k+1 : k-1);
  k=k-1;otherh=otherh-1; // move to c++ indexing
  List ldw_w=as<List>(dw["w"]), ltw_w=as<List>(tw["w"]);
  IntegerVector dw_w(NUMP); // many g will be blank, if obs then length is NUMP
  IntegerVector tw_w(NUMA); // many g will be blank, if obs then length is NUMA
  IntegerVector umatch(max_umatch_size); // many g will be blank, if obs then length is NUMA
  IntegerVector dw_size=as<IntegerVector>(dw["size"]);
  // flip phase at gl if flip[[]][gl]
  for (lk=0;lk<NUMP;lk++) 
    label[lk]=label[lk]-1;
  if (gl==0) {
    h=(flips[gl] ? otherh : k); // if flips[gl]==0 then h=same as was else with even k to k-1 and odd k to k+1
    dw_w=as<IntegerVector>(ldw_w[gl]);
    tw_w=as<IntegerVector>(ltw_w[gl]);
    umatch=as<IntegerVector>(UMATCH[gl]);
    ndonors=kndonors[gl];
    gm=gl*maxdonors*L;
    for (dlk=0;dlk<ndonors;dlk++) 
    {
      dg=gl*offset+dlk; 
      lk=donates[dg];
      llk=label[lk];
      for (l=0;l<L;l++) {
	lm=l*maxdonors;
	fors[gm+lm+dlk]=initProb[l*kLL+llk];
        if (gobs[gl]>0) 
        {
	  gum=(tw_w[h])*dw_size[gl]+dw_w[lk];
	  glk=umatch[gum]; 
	  fors[gm+lm+dlk]*=mutmat[maxmiss*L*glk+L*(gobs[gl]-glk)+l];
        }
	tmpsum+=fors[gm+lm+dlk];
      }
    }
    scalefactor[gl]=1.0/tmpsum;
    for (l=0;l<L;l++) {
      sumfors[gl*L+l]=0.0;
      lm=l*maxdonors;
      for (dlk=0;dlk<ndonors;dlk++) {
	gld=gm+lm+dlk;
	fors[gld]*=scalefactor[gl];
	sumfors[gl*L+l]+=fors[gld];
      }
    }
  }
  if (gl==0) gl=1; // have done initial point above so go to next; else start here
  //for (i=0;i<L;i++) printf("%0.3lf\n", sumfors[1]);
  for (g=gl;g<gu;g++) {
    dw_w=as<IntegerVector>(ldw_w[g]);
    tw_w=as<IntegerVector>(ltw_w[g]);
    umatch=as<IntegerVector>(UMATCH[g]);
    ndonors=kndonors[g-1];
    gm=g*maxdonors*L;
    h=(flips[g] ? otherh : k); // if flips[g]==0 then h=same as was else with even k to k-1 and odd k to k+1
    ndonors=kndonors[g];
    tmpsum=0.0;
    for (dlk=0;dlk<ndonors;dlk++) {
      dg=g*offset+dlk; 
      lk=donates[dg]; 
      llk=label[lk];
      lk2=donatesl[dg]; 
      for (l=0;l<L;l++) {
	lm=l*maxdonors;
	tmpsum2=0.0; // one per l and lk
	lkl=2*L*(llk*L + l); // transition index
	for (i=0;i<L;i++) {
	  ilk=lkl + i; // transition index
	  tmpsum2+=sumfors[(g-1)*L+i]*transitions[ilk];// add in all switch terms
	  if (lk2!=-1) // no match indicated with a -1
	    tmpsum2+=fors[gm-maxdonors*L+i*maxdonors+lk2]*(transitions[ilk+L]-transitions[ilk]); // (switch-to-self+no-switch) and remove already added switch-to-self
	}
	gld=gm+lm+dlk;
	fors[gld]=tmpsum2; // l appears here in case we want anc-specific theta
        if (gobs[g]>0) 
        {
	  gum=(tw_w[h])*dw_size[g]+dw_w[lk];
  	  glk=umatch[gum]; 
	  fors[gld]*=mutmat[maxmiss*L*glk+L*(gobs[g]-glk)+l]; // l appears here in case we want anc-specific theta
        }
	//if (fors[gld]<1.0e-16) fors[gld]=1.0e-16;
	tmpsum+=fors[gld];
      }
    }
    scalefactor[g]=1.0/tmpsum;
    for (l=0;l<L;l++) {
      sumfors[g*L+l]=0.0;
      lm=l*maxdonors;
      for (dlk=0;dlk<ndonors;dlk++) {
	gld=gm+lm+dlk;
	fors[gld]*=scalefactor[g];
	sumfors[g*L+l]+=fors[gld];
      }
    }
  }
}
