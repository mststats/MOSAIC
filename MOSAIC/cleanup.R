# function to clean up files created by MOSAIC that are no longer needed
cleanup_ff_files=function(t.donates, t.donatesl, t.donatesr, t.nchrno, t.NUMI, t.ffpath, t.prethin=FALSE, 
			  t.prethin_donates=NULL,t.prethin_donatesl=NULL,t.prethin_donatesr=NULL)
{
  cat("removing all ff files stored in", t.ffpath, "\n")
  for (ch in 1:t.nchrno) 
    for(ind in 1:t.NUMI) 
    {
      delete(t.donates[[ch]][[ind]])
      delete(t.donatesl[[ch]][[ind]])
      delete(t.donatesr[[ch]][[ind]])
      if (t.prethin)
      {
	delete(t.prethin_donates[[ch]][[ind]])
	delete(t.prethin_donatesl[[ch]][[ind]])
	delete(t.prethin_donatesr[[ch]][[ind]])
      }
    }
}
