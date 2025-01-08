# Polytect

## Introduction
Digital PCR (dPCR) is a state-of-the-art quantification method of target nucleic acids. The technology is based on massive partitioning of a reaction mixture into individual PCR reactions. This results in partition-level end-point fluorescence intensities that are subsequently used to classify partitions as positive or negative, i.e., containing or not containing the target nucleic acid(s). Many automatic dPCR partition classification methods have been proposed, but they are limited to the analysis of single or dual color dPCR data. While general-purpose or flow cytometry clustering methods can be directly applied to multi-color dPCR data, these methods do not exploit the approximate prior knowledge on cluster center locations available in dPCR data.

We developed a novel method, `Polytect`, that relies on crude cluster results from flowPeaks, previously shown to offer good partition classification performance, and subsequently refines flowPeaks' results by cluster merging and automatic cluster labeling, exploiting the prior knowledge on cluster center locations.

The following instructions are for the R package. You can also find the shiny app here: https://digpcr.shinyapps.io/Polytect/

## Installation
`Polytect` can be downloaded from Github.

```r
#GitHub installation
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("emmachenlingo/Polytect")
```

## Examples
To demonstrate the core functions of this package, we will use the HR digital PCR dataset, which can be accessed using the command data(HR). The HR dataset is a data frame consisting of 18,233 rows and 2 columns. The dataset has 4 clusters. The clustering analysis can be performed using the following commands.

```{r load package and data,warning=FALSE, message=FALSE}
library(Polytect)
library(ggplot2)
library(flowPeaks)
data(HR)
head(HR)
ggplot(data=HR, aes(channel1, channel2))+ geom_point(size=0.9,show.legend = FALSE) +
  labs(x = "color 1", y = "color 2") +theme(text = element_text(size = 15), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),
                                                  axis.line = element_line(colour = "black"),
                                                  plot.margin = margin(t = 20,  # Top margin
                                                                       r = 20,  # Right margin
                                                                       b = 20,  # Bottom margin
                                                                       l = 20),
                                                  axis.title.x = element_text(margin = margin(t = 20)),
                                                  axis.title.y = element_text(margin = margin(r = 20)))

```
<div style="text-align: center;">
  <img src="inst/demo_data.png" width="50%" height="70%">
</div>

If we perform flowPeaks only, there will be more clusters than expected.
```{r flowpeaks,warning=FALSE, message=FALSE}
data_scaled<-apply(HR,2,function(x) (x-min(x))/(max(x)-min(x)))
data_input<-as.matrix(data_scaled)
fp<-flowPeaks(data_input)

ggplot(data=HR, aes(channel1, channel2,colour = factor(fp$peaks.cluster)))+ geom_point(size=0.9,show.legend = FALSE) +labs(x = "color 1", y = "color 2") +theme(text = element_text(size = 15), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),
                                                  axis.line = element_line(colour = "black"),
                                                  plot.margin = margin(t = 20,  # Top margin
                                                                       r = 20,  # Right margin
                                                                       b = 20,  # Bottom margin
                                                                       l = 20),
                                                  axis.title.x = element_text(margin = margin(t = 20)),
                                                  axis.title.y = element_text(margin = margin(r = 20)))
```
<div style="text-align: center;">
  <img src="inst/demo_data_fp.png" width="50%" height="70%">
</div>

The main function that performs flowPeaks first, then merges the excess clusters.
```{r flowpeaks and merge,warning=FALSE, message=FALSE}
result<-polytect_clust(data=HR,cluster_num=4)
print(head(result))

ggplot(data=HR, aes(channel1, channel2,colour = factor(result$cluster)))+ geom_point(size=0.9,show.legend = FALSE) +labs(x = "color 1", y = "color 2") +theme(text = element_text(size = 15), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),
                                                  axis.line = element_line(colour = "black"),
                                                  plot.margin = margin(t = 20,  # Top margin
                                                                       r = 20,  # Right margin
                                                                       b = 20,  # Bottom margin
                                                                       l = 20),
                                                  axis.title.x = element_text(margin = margin(t = 20)),
                                                  axis.title.y = element_text(margin = margin(r = 20)))
```
<div style="text-align: center;">
  <img src="inst/demo_data_polytect.png" width="50%" height="70%">
</div>

Or you can use any initial clustering results as an input to the polytect_merge function
```{r merge the data,warning=FALSE, message=FALSE}
## it is advised to standardize the data
dist_matrix <- dist(data_input)
hc <- hclust(dist_matrix, method = "ward.D2")
# the number of clusters is specified at 6, which is larger than 4
hc_clusters <- cutree(hc, k = 6)

ggplot(data=HR, aes(channel1, channel2,colour = factor(hc_clusters)))+ geom_point(size=0.9,show.legend = FALSE) +labs(x = "color 1", y = "color 2") +theme(text = element_text(size = 15), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),
                                                  axis.line = element_line(colour = "black"),
                                                  plot.margin = margin(t = 20,  # Top margin
                                                                       r = 20,  # Right margin
                                                                       b = 20,  # Bottom margin
                                                                       l = 20),
                                                  axis.title.x = element_text(margin = margin(t = 20)),
                                                  axis.title.y = element_text(margin = margin(r = 20)))

hc_parse<-list()
hc_parse$cluster<-hc_clusters

result<-polytect_merge(data=HR,cluster_num=4,base_clust=hc_parse)
print(head(result))
```
<div style="text-align: center;">
  <img src="inst/demo_data_hc.png" width="50%" height="70%">
</div>

The clustering results can be visualized by 2-d plots.
```{r plot the data,warning=FALSE, message=FALSE}
polytect_plot(result,cluster_num=4)
```

<div style="text-align: center;">
  <img src="inst/demo_data_hc_merge.png" width="50%" height="70%">
</div>

You can also summarise the results, which will give you cluster centers, group sizes and silhouette coefficients for each group
```{r summarise the results,warning=FALSE, message=FALSE}
result_summary<-polytect_summary(result)
print(result_summary)
```
<div style="text-align: center;">
  <img src="inst/demo_data_hc_merge.png" width="50%" height="70%">
</div>

You can also plot the individual silhouette coefficient in each cluster
```{r plot sil coefs,warning=FALSE, message=FALSE}
sil_plot(result)
```
<div style="text-align: center;">
  <img src="inst/demo_data_sil.png" width="50%" height="70%">
</div>

There is a function to calculate the concentration of the targets.
```{r calculate the conc, warning=FALSE, message=FALSE}
target_conc<-conc_cal(result,sampvol=0.91,volmix=20,voltemp=20,type="2color")
print(target_conc)
```
<div style="text-align: center;">
  <img src="inst/demo_data_conc.png" width="50%" height="70%">
</div>

This package can also handle 3-up to 6-color dPCR data. We first perform flowPeaks only.
```{r 3d example, warning=FALSE, message=FALSE}
data(BPV)
data_scaled<-apply(BPV,2,function(x) (x-min(x))/(max(x)-min(x)))
data_input<-as.matrix(data_scaled)
fp<-flowPeaks(data_input)
table(fp$peaks.cluster)
df_data<-as.data.frame(cbind(BPV,cluster=fp$peaks.cluster))
polytect_plot(df_data,cluster_num=8)
```
<div style="text-align: center;">
  <img src="inst/bpv_fp_cluster.png" width="50%" height="70%">
</div>

<div style="text-align: center;">
  <img src="inst/bpv_fp.png" width="50%" height="70%">
</div>

Then the main function.
```{r 3d example polytect clust, warning=FALSE, message=FALSE}
result<-polytect_clust(data=BPV,cluster_num=8)
table(result$cluster)
polytect_plot(result,cluster_num=8)
```
<div style="text-align: center;">
  <img src="inst/bpv_polytect_cluster.png" width="50%" height="70%">
</div>

<div style="text-align: center;">
  <img src="inst/bpv_polytect.png" width="50%" height="70%">
</div>



## References

Ge Y, Sealfon S (2012). “flowPeaks: a fast unsupervised clustering for flow cytometry data via K-means and density peak finding.” Bioinformatics. R package version 4.4.0.

Lun A (2024). bluster: Clustering Algorithms for Bioconductor. R package version 1.14.0.
