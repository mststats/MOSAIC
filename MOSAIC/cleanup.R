if (ffcleanup & HPC)
{
  cat("removing all ff files stored in", ffpath, "\n")
  for (ch in 1:nchrno) 
    for(ind in 1:NUMI) 
    {
      delete(donates[[ch]][[ind]])
      delete(donatesl[[ch]][[ind]])
      delete(donatesr[[ch]][[ind]])
      if (prethin)
      {
        delete(prethin_donates[[ch]][[ind]])
        delete(prethin_donatesl[[ch]][[ind]])
        delete(prethin_donatesr[[ch]][[ind]])
      }
    }
}
