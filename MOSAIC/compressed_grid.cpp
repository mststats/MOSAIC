// function to compress genomes to an evenly (on genetic distance) grid. Lists of unique haplotypes at each gridpoint are
// returned along with a local map from each genome to the unique haplotype list. Stored such that re-phasing can be done without recomputing the compression
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List cpp_unique_haps(IntegerVector Y, int k, int S, int G, IntegerVector gmap, int maxgmap, List uhaps) // S is ncol(Y)
{ 
  int s,g,nobs=0,miss=0,l=0,i;
  int mins=0; // loop from mins to S for each g and update as required
  IntegerVector gm(maxgmap); 
  List u; // the unique haps at a gridpoint
  List w; // which haps are identical to which unique haps at each gridpoint
  List du;// decimal version of haplotypes at gridpoints; dy=sum_(i=1)^nobs y_i*2^i
  int gdu=0;
  IntegerVector size; // number of unique haps at this gridpoint
  if (k==1)
  {
    u=List(G); // the unique haps at a gridpoint
    w=List(G); // which haps are identical to which unique haps at each gridpoint
    du=List(G);// decimal version of haplotypes at gridpoints; dy=sum_(i=1)^nobs y_i*2^i
    size=IntegerVector(G);
    for (g=0;g<G;g++) 
    {
      for (s=0;s<maxgmap;s++)
	gm[s]=-1; // no loci at each gridpoint until we find some
      nobs=0;
      u[g]=List(1);
      size(g)=1; // 0=1, etc
      w[g]=IntegerVector(1);as<IntegerVector>(w[g])[0]=0; // 
      for (s=mins;s<S;s++) // loop over possible loci mapped to this gridpoint
      {
	if (gmap(s)>g)
	{
	  mins=s;
	  break;
	}
	if (gmap(s)==g) // this locus maps to this gridpoint
	{
	  gm[nobs]=s; // add to vector here
	  nobs++;
	}
      }
      as<List>(u[g])[0]=IntegerVector(nobs);
      du[g]=IntegerVector(1);as<IntegerVector>(du[g])[0]=0;
      for (s=0;s<nobs;s++)
      {
	as<IntegerVector>(as<List>(u[g])[0])[s]=Y(gm(s)); // add to existing unique haps at this gridpoint
	as<IntegerVector>(du[g])[0]+=(Y(gm(s))*pow(2,s));
      }
    }
  }
  if (k>1)
  {
    u=uhaps["u"]; // already a list of length G
    w=uhaps["w"]; // already a list of length G
    du=uhaps["du"];
    size=as<IntegerVector>(uhaps["size"]);
    for (g=0;g<G;g++) 
    {
      miss=1;
      nobs=0;
      for (s=mins;s<S;s++) // loop over possible loci mapped to this gridpoint
      {
	if (gmap(s)>g)
	{
	  mins=s;
	  break;
	}
	if (gmap(s)==g) // this locus maps to this gridpoint
	{
	  gm[nobs]=s; // add to vector here
	  nobs++;
	}
      }
      if (nobs>0)
      {
	for (l=0;l<size[g];l++)
	{
	  gdu=0;
	  for (s=0;s<nobs;s++)
	    gdu+=(Y(gm(s))*pow(2,s));
	  if (as<IntegerVector>(du[g])[l]==gdu)
	  {
	    miss=0;
	    break; // found a match in previous unique vectors
	  }
	}
	{ // scope so that tmp_w and tmp_ws remove automatically
	  IntegerVector tmp_w(k);
	  IntegerVector tmp_ws=as<IntegerVector>(w[g]);
	  for (i=0;i<(k-1);i++) tmp_w[i]=tmp_ws[i];
	  w[g]=tmp_w;
	}
	if (miss==0) as<IntegerVector>(w[g])[k-1]=l; // seen at l
	if (miss!=0)  // unseen before
	{
	  size[g]++;
	  { // scope so that tmp_u, tmpdu, tmp_us, tmp_dus all remove automatically
	    List tmp_u(size[g]);
	    List tmp_us=as<List>(u[g]);
	    IntegerVector tmp_du(size[g]);
	    IntegerVector tmp_dus=as<IntegerVector>(du[g]);
	    for(i=0; i<size[g]-1; i++) 
	    {
	      tmp_u[i]=tmp_us[i];
	      tmp_du[i]=tmp_dus[i];
	    }
	    IntegerVector newobs(nobs); 
	    for (s=0;s<nobs;s++)
	    {
	      newobs[s]=Y(gm(s)); // add to existing unique haps at this gridpoint
	      tmp_du[size[g]-1]=gdu;
	    }
	    tmp_u[size[g]-1]=newobs; 
	    u[g]=tmp_u;
	    du[g]=tmp_du;
	    as<IntegerVector>(w[g])[k-1]=size[g]-1;
	  }
	}
      }
    }
  }
  List ret;
  ret["w"] = w;
  ret["u"] = u;
  ret["du"] = du;
  ret["size"] = size;
  return wrap(ret);
}
// don't forget to -1 from gmap when calling this...

