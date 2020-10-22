# NPVAR: Learning nonparametric causal graphs in polynomial-time

This is an `R` implementation of the following paper:

[1] Gao, M., Ding, Y., Aragam, B. (2020). [A polynomial-time algorithm for learning nonparametric causal graphs](https://arxiv.org/abs/2006.11970) ([NeurIPS 2020](https://nips.cc/Conferences/2020/)).

If you find this code useful, please consider citing:
```
@inproceedings{gao2020npvar,
    author = {Gao, Ming and Ding, Yi and Aragam, Bryon},
    booktitle = {Advances in Neural Information Processing Systems},
    title = {{A polynomial-time algorithm for learning nonparametric causal graphs}},
    year = {2020}
}
```

## Introduction

The `NPVAR` algorithm is an algorithm for learning the structure of a directed acyclic graph (DAG) `G` that represents a potentially high-dimensional joint distribution `P(X)`. In other words, given `n` samples from `P(X)`, `NPVAR` learns the DAG `G` that generated these samples. In general, this problem is not well-defined since the DAG is not unique, however, we consider the setting where the residual variances `E[var(Xj|pa(j))]` are all approximately equal. In this setting, `NPVAR` is a polynomial-time algorithm for provably recovering the DAG `G`.

## Requirements
- R
- Package `np`
- Package `mgcv`
- Package `igraph`

## Contents
- `NPVAR.R` Main function to run our algorithm, see demo below
- `utils.R` Some helper functions to simulate data and evaluate results
- `ANM_gp.R` File used to generate data from a Gaussian process model, see references

## Demo
Generate a *ER* graph with 5 nodes and 5 expected edges. Then generate data by a *SIN* model with noise variance *sigma*=0.5 and sample size 1000.
```r
source('NPVAR.R')
source('utils.R')

data = data_simu(graph_type = 'ER-SIN', errvar = 0.5, d = 5, n = 1000, s0 = 1, x2 = T)
X = data$X
G = data$G
X2 = data$X2
```
Apply `NPVAR` through 3 implementations:
- Naively recover ordering node by node
- Recover node layer by layer with fixed `eta`
- Recover node layer by layer using `X2` to determine `eta` adaptively
```r
result1 = NPVAR(X)
result2 = NPVAR(X, layer.select = T, eta =  0.01)
result3 = NPVAR(X, layer.select = T, x2 = X2)
```
Check outputs
```r
print(result1)

print(result2$ancestors)
print(result2$layers)

print(result3$ancestors)
print(result3$layers)
```
Furthermore, infer adjacency matrix by significance given by `gam` 
```r
est = prune(X, result1)
print(est)
print(sum(abs(est - G)))
```

## References
- We generate data from Gaussian process models through the `RESIT` code, from [here](https://staff.fnwi.uva.nl/j.m.mooij/code/codeANM.zip)
- Original equal variance code for linear models is from [here](https://github.com/WY-Chen/EqVarDAG)
