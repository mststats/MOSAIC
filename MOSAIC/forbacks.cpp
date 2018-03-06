// this performs the forward-backward algorithm on the MOSAIC HMM (see gfbs)
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector cppforback(int maxdonors,bool THIN,int NUMP,int L,int G,IntegerVector ndonors,IntegerVector donates,NumericVector fors,NumericVector backs)
{
  int offset=(THIN ? maxdonors : 0); // setting offset=0 for maxdonors==NUMP => use single vector of donors
  int NNL2=maxdonors*L, g, l, n, ln, dln;
  double invsum;
  IntegerVector dim(2);
  dim[0]=G; dim[1]=NUMP*L;
  Dimension d(dim);
  NumericVector kforbacks(d);
  for (g=0; g<G; g++) {
    invsum=0.0;
    for (l=0; l<L; l++)
      for (n=0; n<ndonors[g]; n++)
      {
	ln=l*NUMP+donates[g*offset+n];
	dln=l*maxdonors+n;
	kforbacks[ln*G+g]=fors[g*NNL2+dln]*backs[g*NNL2+dln];
	invsum+=kforbacks[ln*G+g];
      }
    invsum=1.0/invsum;
    for (l=0; l<L; l++)
      for (n=0; n<ndonors[g]; n++)
      {
	ln=l*NUMP+donates[g*offset+n];
	kforbacks[ln*G+g]*=invsum;
      }
  }
  return(wrap(kforbacks));
}

// [[Rcpp::export]]
NumericVector cppgforback(int maxdonors,int THIN,int kLL,int NUMP,IntegerVector label,int L,int G,IntegerVector ndonors,IntegerVector donates,
    NumericVector fors,NumericVector backs)
{
  int offset=(THIN ? maxdonors : 0); // setting offset=0 for maxdonors==NUMP => use single vector of donors
  int NNL2=maxdonors*L, g, l, n, ln, dln, lk;
  double tmp, invsum;
  IntegerVector dim(2);
  dim[0]=G; dim[1]=kLL*L;
  Dimension d(dim);
  NumericVector kforbacks(d);
  for (lk=0;lk<NUMP;lk++) 
    label[lk]=label[lk]-1;
  for (g=0; g<G; g++) {
    invsum=0.0;
    for (l=0; l<L; l++)
      for (n=0; n<ndonors[g]; n++)
      {
	ln=l*kLL+label[donates[g*offset+n]];
	dln=l*maxdonors+n;
	tmp=fors[g*NNL2+dln]*backs[g*NNL2+dln];
	kforbacks[ln*G+g]+=tmp;
	invsum+=tmp;
      }
    invsum=1.0/invsum;
    for (l=0; l<L; l++)
      for (n=0; n<kLL; n++)
      {
	ln=l*kLL+n;
	kforbacks[ln*G+g]*=invsum;
      }
  }
  return(wrap(kforbacks));
}

