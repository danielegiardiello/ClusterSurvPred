

# The impact of the size and the number of clusters on prediction performance of the stratified and the conditional shared gamma frailty Cox proportional hazards models

**Authors**: Daniele Giardiello, Edoardo Ratti, Peter C. Austin

## Abstract
Researchers in biomedical research often analyse data that are subject to clustering.
Independence among observations are generally assumed to develop and validate risk
prediction models. For survival outcomes, the Cox proportional hazards regression model
is commonly used to estimate individual’s risk at fixed time horizons. The stratified Cox
proportional hazards and the shared gamma frailty Cox proportional hazards regression
models are two common approaches to account for the presence of clustering in the data.
The accuracy of the predictions of these two approaches has not been examined. We
conducted a set of Monte Carlo simulations to assess the impact of the cluster size,
number of clusters, and the within-cluster correlation on the accuracy of the conditional
predictions developed using the stratified and the shared gamma frailty Cox proportional
hazards regression model. We compared the accuracy of the predictions in terms of
discrimination, calibration and overall performance metrics. We found that the stratified
and the shared gamma frailty model provided similar performances especially for larger
size and higher number of clusters. For small cluster size, simulations suggested slightly
better discrimination and overall performance for the stratified model and slightly better
calibration for the shared gamma frailty model at short time horizons. However, the
practical applicability of the stratified Cox proportional hazards model to estimate
predictions is limited especially for high within-cluster correlation and when clusters are
small, and more likely at longer time prediction horizons.


## Usage
The illustrative case study accompanying this work is the insem data from the R package [`parfm`](https://cran.r-project.org/web/packages/parfm/index.html)

You can either download a zip file containing the directory, or you can clone it by using

```bash
git clone https://github.com/danielegiardiello/ClusterSurvPred.git
```

In either case, you can then use the `ClusterSurvPred.Rproj` file to open
and Rstudio session in the directory you have just downloaded. You may then knit
both rmarkdown files, or run them line-by-line.

## Contributions

| Name                                                         | Affiliation                           | Role                  |
| ------------------------------------------------------------ | ------------------------------------- | ----------------------|
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | University of Milan-Bicocca (IT) | Author/maintainer |
| [Edoardo Ratti](https://en.unimib.it/edoardo-ratti) | University of Milan-Bicocca (IT) | Code reviewer        |
| [Peter C. Austin](https://www.ices.on.ca/ices-scientists/peter-austin/) | ICES Toronto (CA)  <br /> Institute of Health Policy, Management and Evaluation, University of Toronto (CA) <br /> Sunnybrook Research Institute, Toronto (CA) | Contributor |





