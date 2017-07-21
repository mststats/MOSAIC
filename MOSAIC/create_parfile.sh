#!/bin/sh
if [ "$#" -ne 10 ]; then
  echo "usage: sh create_parfile.sh alpha lambda rho1 rho2 theta1 theta2 theta3 ww[1,2] ww[2,1] chrno"
  exit
fi 
# alpha=1 lambda=2 rho1=3 rho2=4 theta1=5 theta2=6 theta3=7 ww[1,2]=8 ww[2,1]=9 chrno=9"
cp mosaic.par.txt mosaic_par.${10}
sed -i s/"THETAVAL"/$1/g mosaic_par.${10}
sed -i s/"LAMBDAVAL"/$2/g mosaic_par.${10}
sed -i s/"RHO1VAL"/$3/g mosaic_par.${10}
sed -i s/"RHO2VAL"/$4/g mosaic_par.${10}
sed -i s/"MUT1VAL"/$5/g mosaic_par.${10}
sed -i s/"MUT2VAL"/$6/g mosaic_par.${10}
sed -i s/"MUT3VAL"/$7/g mosaic_par.${10}
sed -i s/"MISCOPY1VAL"/$8/g mosaic_par.${10}
sed -i s/"MISCOPY2VAL"/$9/g mosaic_par.${10}
sed -i s/"CHRVAL"/${10}/g mosaic_par.${10}
# create a version for unphased data also
cp mosaic_par.${10} mosaic_gen_par.${10}
sed -i s/"GENOTYPE:0"/"GENOTYPE:1"/g mosaic_gen_par.${10}
sed -i s/"HAPMIX_MODE:HAPLOID"/"HAPMIX_MODE:LOCAL_ANC"/g mosaic_gen_par.${10}
sed -i s/"admixedhaplofile"/"admixedgenofile"/g mosaic_gen_par.${10}

