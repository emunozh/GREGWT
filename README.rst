========
R GREGWT
========

:Author: Esteban Munoz
:Version: 0.7.0
:Date: Mon Aug 19. 2015

.. contents:: Table of Contents
   :depth: 2

Implementation of GREGWT in R. 

Install the development version on GREGWT
-----------------------------------------

Using the devtools library::

    library(devtools) 
    install_github("emunozh/GREGWT")

Prepare data for simulation
---------------------------

::

    Simulation.Data <- prepareData(
        GREGWT.census, GREGWT.survey,
        # census parameters
        pop_benchmark=c(2,12),
        census_categories=seq(2,24),
        # survey parameters
        survey_id=FALSE,
        survey_categories=seq(1,3)
    )

Runin GREGWT
------------

Rewheight the survey to area ``area_code``::

    Weights.GREGWT <- GREGWT(data_in=Simulation.Data, area_code="02")

Create a synthetic population (under development)
-------------------------------------------------

Prepare a survey to be used by an agent based simulation model::

    Synthetic.pop <- Synthetize(data_in=Weights.GREGWT)


Data:
-----

Two datasets are provided in this package. 

1. Census data from 2011 aggregated to municipalities. (GREGWT.census)::

       data("GREGWT.census")

2. A reformatted version of the scientific data file from the micro-census 2010
   (GREGWT.survey)::

       data("GREGWT.survey")
