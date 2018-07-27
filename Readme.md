# PeriFlakes

Alejandro Francisco Queiruga  
Lawrence Berkeley National Lab  
2015-2018

## Introduction

PeriFlakes is an open source solver of state-based Peridynamics.
It is built off of cornflakes/popcorn and can autogenerate every possible
formulation found in the literature, and then some more. Notably, the system
uses symbolic differentiation to obtain the expressions for tangent matrices to
use implicit time steppers or solve static equilibrium. This is an updated
version of the code used to perform the study in

> Queiruga, A. F. and G. J. Moridis, "Numerical experiments on
>  the convergence properties of state-based peridynamic laws and
>  influence functions in two-dimensional problems." Computer
>  Methods in Applied Mechanics and Engineering 322 (2017):
>  97-122.

in which the algorithms and formulations are described in detail.

## Running

The easiest way to run this code yourself is to pull the cornflakes image:
```python
docker pull afqu/cornflakes
```
and then find this code in `/opt_cornflakes/PeriFlakes`.

## Datasets

The databases that the scripts will generate are published on Zenodo
at:

> Alejandro Francisco Queiruga. (2018). 2D Peridynamic Displacement
> Fields [Data set]. Zenodo. http://doi.org/10.5281/zenodo.1284634

They are simple sqlite3 databases, and the Jupyter notebook
[analysis.ipynb](analysis.ipynb) illustrates their structure and how to query
them.

## Conclusions

The purpose of this repository is to illustrate the limitations of the
Peridynamics family of numerical methods.
The popcorn file
[PeriFlakes/peri_kernels_pop.py](PeriFlakes/peri_kernels_pop.py)
performs extensive Ahead-Of-Time code generation to allow us to
exhaustively search the hyperparameter space of possible Peridynamics programs.
The simulation object was only designed to quickly graph the unit square
domain and sweep through different configurations of the Peridynamics
approximation.

The analysis of the results can be viewed and run in the Jupyter
notebook [analysis.ipynb](analysis.ipynb).
The fundamental limitation is the inability of the most widely used state-based
to represent all basic modes of deformation of mechanics.
The derfomation-gradient based correspondence model is able to
converge linearly to all explored problems with certain hyperparameter choices,
and its instability can be  dealt with.
The recommendation for using perdiynamics is thus:

- Use the deformation gradient based correspondence model of Silling,2007 using
either the stabilizer of Silling,2017 or the post-smoothing technique of
Queiruga, 2017.
- Use a cubic influence function with as small of an influence radius as
possible, no shorter than $1.5\times$ the grid size.

Some arguments on the computational complexity of Peridynamics are made in the
Jupyter notebook [computation.ipynb](computation.ipynb).

The ficticious-nodes boundary correction method is explored in
[ficticious_nodes.ipynb](ficticious_nodes.ipynb).

## Acknowledgements

This work was supported primarily by Laboratory Directed Research and
Development (LDRD) funding from Berkeley Lab, provided by the Director, Office
of Science, of the U.S. Department of Energy under Contract No.
DE-AC02-05CH11231.

## License

Copyright (C) Alejandro Francisco Queiruga, LBNL, 2015-2018

PeriFlakes is released under version 3 of the GNU General Public License, as per
LICENSE.txt.
