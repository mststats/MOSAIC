#!/bin/sh
if [ "$#" -ne 1 ]; then
  echo "Usage: location_to_install"
  exit 1
fi
rlib=$1
rm -rf MOSAIC/ MOSAIC_1.0.tar.gz 
Rscript create_package.R
echo "importFrom(parallel,detectCores)" >> MOSAIC/NAMESPACE
echo "importFrom(doParallel,registerDoParallel)" >> MOSAIC/NAMESPACE
echo "importFrom(foreach,foreach)" >> MOSAIC/NAMESPACE
echo "importFrom(foreach,'%dopar%')" >> MOSAIC/NAMESPACE
echo "importFrom(ff,ff)" >> MOSAIC/NAMESPACE
echo "importFrom(ff,delete)" >> MOSAIC/NAMESPACE
echo "importFrom(compiler,cmpfun)" >> MOSAIC/NAMESPACE
echo "importFrom(cluster,fanny)" >> MOSAIC/NAMESPACE
echo "importFrom(combinat,permn)" >> MOSAIC/NAMESPACE
#echo "importFrom(argparser,arg_parser)" >> MOSAIC/NAMESPACE
#echo "importFrom(argparser,add_argument)" >> MOSAIC/NAMESPACE
#echo "importFrom(argparser,parse_args)" >> MOSAIC/NAMESPACE
cp DESCRIPTION MOSAIC/DESCRIPTION

cp -r README.txt mosaic.R example_data/ MOSAIC/ 

rm MOSAIC/Read-and-delete-me

R CMD build MOSAIC
R CMD INSTALL -l $rlib MOSAIC_1.0.tar.gz

