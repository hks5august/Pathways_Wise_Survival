c=".csv"
while read line
do
mkdir $line

cd $line
wait 

cp ../Univariate_survival_analysis_script.R .
cp ../$line$c .

wait
nohup Rscript Univariate_survival_analysis_script.R $line$c &

wait
cd ../

done<list

