while read line
do

cd $line
cp ../*.R .

Rscript multivariate_with_PI_script.R 

cd ../

done <list_10
