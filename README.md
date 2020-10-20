# NPVAR
 NPVAR: Learning Nonparametric DAG with variance assumption

This is an implementation of paper: 
[A polynomial-time algorithm for learning nonparametric causal graphs](https://arxiv.org/abs/2006.11970)

## Requiment
- R
- Package `np`
- Package `mgcv`
- Package `igraph`

## Contents
- `NPVAR.R` Main function to run our algorithm, see demo below
- `utils.R` Some helper function to simulate data and evaluate results
- `ANM_gp.R` File used to generate Gaussian process model, see reference

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

## Reference
We generate Gaussian process model through code of `RESIT` from [here](https://staff.fnwi.uva.nl/j.m.mooij/code/codeANM.zip)

Naive equal variance code is from [here](https://github.com/WY-Chen/EqVarDAG)
