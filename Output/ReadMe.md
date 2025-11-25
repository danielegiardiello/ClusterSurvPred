# Output folder

## Proposed case study

+ df_perf_casestudy_cimbtr.rds provides the accuracy of the prediction for the CIMBTR data case study about mortality after transplant among patients with myelodysplastic syndrome.  

+ df_perf_casestudy_bladder.rds provides the accuracy of the prediction for the bladder cancer data.  



## Extra case study insem data

+ df_perf_casestudy_large.rds provides the accuracy of the predictions for the large cluster scenario of the illustrative case study 
over the 100 repeated random splits (70% derivation, 30% validation) stratified by the cluster

+ df_perf_casestudy_small.rds provides the accuracy of the predictions for the small cluster scenario of the illustrative case study 
over the 100 repeated random splits (70% derivation, 30% validation) stratified by the cluster


The large cluster scenario included at least 50 subjects per cluster (i.e., 50 cows per herd).  

The small cluster scenario included less than 50 subjects per cluster (i.e, less than 50 cows per herd).  
Herds with less than 10 cows per herd were excluded due to the difficulties to obtain the same herds in both the derivation and validation samples
