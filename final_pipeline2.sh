#a="pathway_sig_genes.txt"
a="New.tsv"
while read line 

do

echo $line
cd $line/
cp ../*.R .

wait

nohup Rscript pathways_prognostic_index_models2.R  $line$a &
#nohup   Rscript selected_feature_data_prep_script.R $line$a sel_train sel_test sel_ext  &
wait

perl ~/scripts/hashmatch.pl tr_PI ../tr_clin 1 1  |cut -f1-25,28 >tr_clin_PI

perl ~/scripts/hashmatch.pl te1_PI ../te_clin 1 1  |cut -f1-25,28 >te1_clin_PI

perl ~/scripts/hashmatch.pl te2_PI ../ext_clin 1 1  |cut -f1-25,28 >te2_clin_PI

wait

Rscript multivariate_with_PI_script.R 

wait

cd ../

done<25_list
