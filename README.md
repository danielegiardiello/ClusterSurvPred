

# The impact of the number and the size of clusters on prediction performance of the stratified and the conditional shared gamma frailty Cox proportional hazards models

[![medRxiv doi](https://img.shields.io/badge/medRxiv–doi-10.1101%2F2025.10.17.25338219v1-FFA500.svg)](https://doi.org/10.1101/2025.10.17.25338219)

**Authors**: Daniele Giardiello, Edoardo Ratti, Peter C. Austin


## Abstract
 <div align="justify"> 
Researchers in biomedical research often analyse data that are subject to clustering. Independence among observations are generally assumed to develop and validate risk prediction models. For survival outcomes, the Cox proportional hazards regression model is commonly used to estimate an individual’s risk at fixed time horizons. The stratified Cox proportional hazards and the shared gamma frailty Cox proportional hazards regression models are two common approaches to account for the presence of clustering in the data. The accuracy of the predictions of these two approaches has not been examined. We conducted a set of Monte Carlo simulations to assess the impact of the number of clusters, the size of the clusters, and the within-cluster correlation in outcomes on the accuracy of the conditional predictions developed using the stratified and the shared gamma frailty Cox proportional hazards regression model. We compared the accuracy of the predictions in terms of discrimination, calibration and overall performance metrics. We found that the stratified and the shared gamma frailty model provided similar performance, especially for larger size and higher number of clusters. For small cluster size, we observed slightly better discrimination and overall performance for the stratified model and better calibration for the shared gamma frailty model at shorter prediction horizons. However, the practical applicability of the stratified Cox proportional hazards model to estimate predictions is limited especially for high within-cluster correlation and when clusters are small, and more likely at longer time prediction horizons.
</div>

## Usage
The illustrative case study accompanying this work is the insem data from the R package [`parfm`](https://cran.r-project.org/web/packages/parfm/index.html)

You can either download a zip file containing the directory, or you can clone it by using

```bash
git clone https://github.com/danielegiardiello/ClusterSurvPred.git
```

In either case, you can then use the `ClusterSurvPred.Rproj` file to open
and Rstudio session in the directory you have just downloaded. You may then knit
both rmarkdown files, or run them line-by-line.

Code for the illustrative case study: 
+ Directly available online using [01_casestudy_md.md](https://github.com/danielegiardiello/ClusterSurvPred/blob/main/01_casestudy_md.md) or using the [pdf version](https://github.com/danielegiardiello/ClusterSurvPred/blob/main/01_casestudy_pdf.pdf)
  
+  One can download and open the [html file](https://github.com/danielegiardiello/ClusterSurvPred/blob/main/01_casestudy_html.html)

+ The corresponding source files are the .Rmd files. We define three different source files since there are some cosmetical changes. 

``` r
.
├── 01_casestudy_html.html
├── 01_casestudy_html.Rmd
├── 01_casestudy_md.md
├── 01_casestudy_md.Rmd
├── 01_casestudy_pdf.pdf
├── 01_casestudy_pdf.Rmd
├── ClusterSurvPred.Rproj
├── Data
│   └── insem.rds
├── Functions
│   ├── NEWS.md
│   ├── predict.coxph.gammafrail.R
│   ├── rate_cens_bisection.R
│   └── ReadMe.md
├── imgs
│   ├── additional-1.pdf
│   ├── additional-1.png
│   ├── descriptives-1.pdf
│   ├── descriptives-1.png
│   ├── descriptives-2.pdf
│   ├── descriptives-2.png
│   ├── KMcurves-1.png
│   └── plot_cluster_size-1.png
├── LICENSE
├── Output
│   ├── df_perf_casestudy_large.rds
│   ├── df_perf_casestudy_small.rds
│   └── ReadMe.md
└── README.md
```

## Data availability statement
The first case study dataset was collected by the Center for International Blood and Marrow Transplant
Research (CIBMTR) which is supported primarily by the Public Health Service U24CA076518 from the
National Cancer Institute; the National Heart, Lung, and Blood Institute; the National Institute of Allergy
and Infectious Diseases; 75R60222C00011 from the Health Resources and Services Administration; N00014-
23-1-2057 and N00014-24-1-2507 from the Office of Naval Research; NMDP; and the Medical College of
Wisconsin.

## Contributions

| Name                                                         | Affiliation                           | Role                  |
| ------------------------------------------------------------ | ------------------------------------- | ----------------------|
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | University of Milan-Bicocca (IT) | Author/maintainer |
| [Edoardo Ratti](https://en.unimib.it/edoardo-ratti) | University of Milan-Bicocca (IT) | Code reviewer        |
| [Peter C. Austin](https://www.ices.on.ca/ices-scientists/peter-austin/) | ICES Toronto (CA)  <br /> Institute of Health Policy, Management and Evaluation, University of Toronto (CA) <br /> Sunnybrook Research Institute, Toronto (CA) | Contributor |





