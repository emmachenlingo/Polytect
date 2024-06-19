# Polytect

## Introduction
Digital PCR (dPCR) is a state-of-the-art quantification method of target nucleic acids. The technology is based on massive partitioning of a reaction mixture into individual PCR reactions. This results in partition-level end-point fluorescence intensities that are subsequently used to classify partitions as positive or negative, i.e., containing or not containing the target nucleic acid(s). Many automatic dPCR partition classification methods have been proposed, but they are limited to the analysis of single or dual color dPCR data. While general-purpose or flow cytometry clustering methods can be directly applied to multi-color dPCR data, these methods do not exploit the approximate prior knowledge on cluster center locations available in dPCR data.

We developed a novel method, `Polytect`, that relies on crude cluster results from flowPeaks, previously shown to offer good partition classification performance, and subsequently refines flowPeaks' results by cluster merging and automatic cluster labeling, exploiting the prior knowledge on cluster center locations.


## Installation
`Polytect` can be directly downloaded from Bioconductor.

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("Polytect")

#GitHub installation
install.packages("devtools")
library(devtools)
install_github("emmachenlingo/Polytect")
```

## Examples
To demonstrate the core functions of this package, we will use the HR digital PCR dataset, which can be accessed using the command data(HR). The HR dataset is a data frame consisting of 18,233 rows and 2 columns. The clustering analysis can be performed using the following commands.

```{r load package and data,warning=FALSE, message=FALSE}
library(Polytect)
data(HR)
head(HR)
plot(HR)
```

The main function that performs flowPeaks first, then merges the excess clusters.
```{r merge the data,warning=FALSE, message=FALSE}
result<-polytect_merge(data=HR,cluster_num=4,type="2color")
print(head(result))
```


The clustering results can be visualized by 2-d plots.
```{r plot the data,warning=FALSE, message=FALSE}
polytect_plot(result)
```
You can also summarise the results, which will give you cluster centers, group sizes and silhouette coefficients for each group
```{r summarise the results,warning=FALSE, message=FALSE}
result_summary<-polytect_summary(result)
print(result_summary)
```

You can also plot the individual silhouette coefficient in each cluster
```{r plot sil coefs,warning=FALSE, message=FALSE}
sil_plot(result)
```

There is a function to calculate the concentration of the targets.
```{r calculate the conc, warning=FALSE, message=FALSE}
target_conc<-conc_cal(result,sampvol=0.91,volmix=20,voltemp=20,type="2color")
print(target_conc)
```

## Session Information

The following information was obtained from the R session that generated this vignette:

```{r session_info, eval = FALSE}
sessionInfo()
```

## References

Ge Y, Sealfon S (2012). “flowPeaks: a fast unsupervised clustering for flow cytometry data via K-means and density peak finding.” Bioinformatics. R package version 4.4.0.

Lun A (2024). bluster: Clustering Algorithms for Bioconductor. R package version 1.14.0.
