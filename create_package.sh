#!/bin/sh
if [ "$#" -ne 1 ]; then
  echo "Usage: location_to_install"
  exit 1
fi
rlib=$1
rm -rf MOSAIC/ MOSAIC_1.0.tar.gz 
Rscript create_package.R
echo "import(parallel)" >> MOSAIC/NAMESPACE
echo "import(compiler)" >> MOSAIC/NAMESPACE
echo "importFrom(cluster, fanny)" >> MOSAIC/NAMESPACE
echo "importFrom(combinat, permn)" >> MOSAIC/NAMESPACE
#echo "import(bit) " >> MOSAIC/NAMESPACE

rm MOSAIC/Read-and-delete-me

R CMD build MOSAIC
R CMD INSTALL -l $rlib MOSAIC_1.0.tar.gz

