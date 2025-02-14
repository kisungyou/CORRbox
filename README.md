# CORRbox

This repository contains **CORRbox**, a MATLAB toolbox that accompanies the manuscript titled "_Scalable Geometric Learning with Correlation-Based Functional Brain Networks_", co-authored with [Prof. Hae-Jeong Park](http://neuroimage.yonsei.ac.kr/). Please note that this toolbox is an enhanced alternative to [another package](https://github.com/kisungyou/papers/tree/master/02-CORRbox) that uses quotient geometry in terms of computational efficiency.

## 1. Requirements and Installation
This package is designed for use with `MATLAB`, requiring some modules such as 

  - [Optimization Toolbox][1]
  - [Statistics and Machine Learning Toolbox][2]

After downloading this repo, place `CORRbox` folder to where you want to use and add it onto your working path in MATLAB.

## 2. How to use?

I added a number of example scripts in the folder with detailed comments. Please consult those. For all functions in this package, the data is assumed to be _wrapped_ using `corr_init` functions, which performs *preparation* of the data by checking whether each FC matrix is truly a correlation matrix.

If you want to find documentation of certain functions, you may type in MATLAB console like typical functions. For instance, suppose you are interested in learning how to use the function `corr_mean`. Then, the help page can be shown by typing the following in the console:

```
help corr_mean
```
If you want a pop-up window for a dedicated help page, type the following in the console:
```
doc corr_mean
```
which will show a separate help page.

## 3. What is geometry?

In many functions, `Name, Value` pairs are available to opt for, where you can find all possible combinations from the help page. In most cases, there will be an option called `geometry`, which takes one of the two values: `"ECM"` or `"LEC"`. They are two supported geometries that will appear in our manuscript. In short, these geometries are why and how algorithms in this package can be boosted in terms of efficiency. 

  
## 4. What kinds of jobs can you do?

There are a number of functions for different tasks. All the functions are named as `corr_*` for consistency. Please see the table below on what it can do.

### (1) low-dimensional embedding/visualization

Given $N$ correlation matrices, find an $N\times p$ matrix where each row represents a low-dimensional embedding of each correlation matrix, with $p=2,3$ for visualization purpose.

| Algorithm | Description |
| :------- | :----------- |
| **`corr_ae`**| autoencoder with a single hidden layer. |
| **`corr_cmds`**| classical multidimensional scaling. |
| **`corr_mmds`**| metric multidimensional scaling. |
| **`corr_pga`**| principal geodesic analysis. |
| **`corr_tsne`**| t-stochastic neighbor embedding. |

### (2) cluster analysis

Given $N$ correlation matrices, cluster analysis aimes to find $K$ clusters. 

| Algorithm | Description |
| :------- | :----------- |
| **`corr_kmeans`**|k-means clustering |
| **`corr_kmedoids`**| k-medoids clustering |
| **`corr_specc`**|spectral clustering |
| **`corr_silhouette`**|cluster validity index of Silhouette score |
| **`corr_CH`**|cluster validity index of Calinski and Harabasz |

### (3) regression analysis

Suppose you have $N$ correlation matrices $X_1, \ldots, X_N$. The class of regression algorithms aims at finding nonlinear relationship to predict scalar-valued responses $Y_1, \ldots, Y_N$.

| Algorithm | Description |
| :------- | :----------- |
| **`corr_gpreg`**| gaussian process regression |
| **`corr_kernreg`**|kernel regression |
| **`corr_svmreg`**|support vector machine regression |


### (4) hypothesis testing

All the routines under this category are testing whether two sets of correlation-valued observations are from same distribution or not. 

| Algorithm | Description |
| :------- | :----------- |
| **`corr_test2mmd`**|two-sample testing via Maximum Mean Discrepancy |
| **`corr_test2bg`**|two-sample testing by Biswas-Ghosh method |
| **`corr_test2wass`**| two-sample testing with Wasserstein distance|
| **`corr_test2energy`**|two-sample testing by Energy distance |

### (5) exploratory analysis

These are supplementary functions that can either be used on their own for exploring the nature of the dataset or basis of other codes.

| Algorithm | Description |
| :------- | :----------- |
| **`corr_mean`**|Fréchet mean |
| **`corr_median`**| Fréchet median|
| **`corr_mmd`**|Maximum Mean Discrepancy |
| **`corr_wasserstein`**|Wasserstein distance between two empirical measures |

  
  
  
  
[1]: https://www.mathworks.com/products/optimization.html
[2]: https://www.mathworks.com/products/statistics.html