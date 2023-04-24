a="pathway_sig_genes.txt"
while read line 

do

echo $line
cd $line/
cp ../*.R .

wait


nohup   Rscript selected_feature_data_prep_script.R $line$a sel_train sel_test sel_ext  &
wait

nohup Rscript All_ML_Model.R sel_train sel_test	Recurrent Train_out_res &
wait

nohup Rscript external_val_new.R sel_test  Recurrent  New_test1_res New_test1_pred New_test1_prob_pred  New_test1_acc  &

wait
nohup Rscript ROC_Value.R New_test1_pred test1_roc &

mv Rplots.pdf test1_ROC_plot.pdf

wait

wait
nohup Rscript external_val_new.R sel_ext  Recurrent  New_test2_res New_test2_pred New_test2_prob_pred  New_test2_acc &

wait

nohup Rscript ROC_Value.R New_test2_pred test2_roc &

 mv Rplots.pdf test2_ROC_plot.pdf

cd ../

done <25_list
