========
R GREGWT
========

:Author: Esteban Munoz
:Version: 1.0
:Date: Mon Oct 20. 2014

.. contents:: Table of Contents
   :depth: 2

Implementation of GREGWT in R. 

How to use:
-----------

See the `./Examples` folder for more information.

Basic example::

    Weights.New = GREGWT(X, dx, Tx, group="HHid", bounds=c(0, Inf))

Where

1. `X` is the sample, formated either as a `matrix` or as a `data.frame`
2. `dx` are the initial weights formated as a vector
3. `Tx` are the true population totals
4. (Optional) `group` can be set to define one of the columns of `X` to set a grouping parameter (e.g. households id's)
5. (Optional) `bounds` sets the truncation bounds as `c(L, U)`. Default values are: `c(-Inf, Inf)`
6. (Optional) `epsilon` defining the convergence criterion. Default is set to: `epsilon = 0.001`.
7. (Optional) `max.iter` defining the maximum number of iterations. Default isset to: `max.iter = 10`.

Data:
-----

Two datasets are provided in this package. 

1. Census data from 2011 aggregated to municipalities. (GREGWT.census)::

       data("GREGWT.census")

2. A reformatted version of the scientific data file from the micro-census 2010
   (GREGWT.survey)::

       data("GREGWT.survey")

The following section briefly describes the data structure of both datasets:
the Census 2011 and the Micro Census Survey 2002.

GREGWT.census
+++++++++++++

GREGWT.survey
+++++++++++++

Prepare data for simulation
---------------------------

::

    prepareData(X, Tx, dx, group=F, reference.col=F, covert=T)




