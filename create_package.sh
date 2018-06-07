#!/bin/sh
if [ "$#" -ne 1 ]; then
  echo "Usage: location_to_install"
  exit 1
fi
rlib=$1
rm -rf mosaicpackage/ mosaicpackage_1.0.tar.gz 
Rscript create_package.R
R CMD build mosaicpackage
R CMD INSTALL -l $rlib mosaicpackage_1.0.tar.gz

