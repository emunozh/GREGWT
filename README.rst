========
R GREGWT
========

:Author: Esteban Munoz
:Version: 0.7.3
:Date: Sat 19 Mar 2016 06:27:35 PM CET

.. contents:: Table of Contents
   :depth: 2

Implementation of GREGWT in R. 

Install the development version on GREGWT
-----------------------------------------

Using the devtools library

.. code:: r

   library(devtools) 
   install_github("emunozh/GREGWT")

Prepare data for simulation
---------------------------

.. code:: r

    Simulation.Data <- prepareData(
        GREGWT.census, GREGWT.survey,
        # census parameters
        pop_benchmark=c(2,12),
        census_categories=seq(2,24),
        # survey parameters
        survey_id=FALSE,
        survey_categories=seq(1,3)
    )

Parameters:
^^^^^^^^^^^

Running GREGWT
--------------

Rewheight the survey to area ``area_code``

.. code:: r

    Weights.GREGWT <- GREGWT(data_in=Simulation.Data, area_code="02")

Create a synthetic population (under development)
-------------------------------------------------

Prepare a survey to be used by an agent based simulation model

.. code:: r

    Synthetic.pop <- Synthetize(data_in=Weights.GREGWT)

Data
----

Two datasets are provided in this package. 

1. Census data from 2011 aggregated to municipalities ``GREGWT.census``.
   Original data set downloaded from [1]_.

   .. code:: r

       data("GREGWT.census")

2. A reformatted version of the public available micro-census of year 2010
   (CAMPUS-File)
   This survey data contains (10179) records and corresponds to the (3.5%)
   of the official micro census data for scientific use (1% of the total
   German population) ``GREGWT.survey``.
   Original data set downloaded from [2]_.

   .. code:: r

       data("GREGWT.survey")

.. [1] https://ergebnisse.zensus2011.de
.. [2] http://www.forschungsdatenzentrum.de/bestand/mikrozensus/cf/2010/fdz_mikrozensus_cf_2010_ascii-csv.zip
