
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ORsurv - Survival of patients using ors on TCGA data.

Author: Krishan Gupta

## Introduction

ORsurv - survival of patients using ors on TCGA data.

## Installation

The developer version of the R package can be installed with the
following R commands:

install\_github(“krishan57gupta/ORsurv”)

## Vignette tutorial

This vignette uses a small simulated data, to demonstrate a standard
pipeline. This vignette can be used as a tutorial as well.

## Example

Libraries need to be loaded before running.

``` r
library(ORsurv)
```

``` r
matrix_1=matrix(-499:500,100,10)
rownames(matrix_1)=paste("R",1:100,sep="")
colnames(matrix_1)=paste("C",1:10,sep="")
matrix_2=matrix(1:1000,100,10)
rownames(matrix_2)=paste("R",1:100,sep="")
colnames(matrix_2)=paste("C",1:10,sep="")
matrix_3=matrix(1:1000,100,10)
rownames(matrix_3)=paste("R",1:100,sep="")
colnames(matrix_3)=paste("C",1:10,sep="")
t1=-3
t2=3
t3=5
method="cosine"
name1=paste("EXP_new",method,t2,t3,sep="_")
name2=paste("EXP_TCGA_new",method,t2,t3,sep="_")
name3=paste("Survival_plot_new",method,t2,t3,sep="_")
survival_info=matrix(c(1:10,rep(1,5),rep(0,5)),10,2)
rownames(survival_info)=paste("C",1:10,sep="")
colnames(survival_info)=paste("R",1:2,sep="")
month_limit=180
p_limit=2
selected_label=c(1,2,3,4,5)
folder="~/Ahuja_Lab/EXP52_3/"
```

``` r
output<-ORsurv(folder=folder,matrix_1=matrix_1,matrix_2=matrix_2,matrix_3=matrix_3,t1=t1,t2=t2,t3=t3,
         name1=name1,name2=name2,name3=name3,survival_info=survival_info,method=method,
         month_limit=month_limit,p_limit=p_limit,selected_label=selected_label)
```

<img src="man/figures/README-main-1.png" width="100%" /><img src="man/figures/README-main-2.png" width="100%" />

``` r
output
#>    strata     median     lower upper
#> 1 label=1 0.03333333        NA    NA
#> 2 label=2 0.06666667        NA    NA
#> 3 label=3 0.10000000        NA    NA
#> 4 label=4 0.16666667 0.1333333    NA
#> 5 label=5         NA        NA    NA
```
