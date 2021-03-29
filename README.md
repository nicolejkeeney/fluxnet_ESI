# Project description

Code, figures, and data for Fall 2020 senior thesis at UC Berkeley in Atmospheric Science, written by Nicole Keeney and advised by Prof Dennis Baldocchi. 

Title: Evaluation of a simple parameterization of the Evaporative Stress Index using FLUXNET data and a planetary boundary layer model

Contact: nicolejkeeney @ gmail.com

# Abstract 

The Evaporative Stress Index (ESI) is a drought index that describes temporal anomalies in evapotranspiration, providing an important measure of ecosystem water stress and soil moisture deficits. We evaluate a simple parameterization of ESI introduced by Fisher et al. (2008) as a soil moisture constraint that is reliant only on relative humidity, vapor pressure deficit, and a power law parameter, greatly reducing the number of input parameters required and allowing for calculation of soil moisture deficits using atmospheric variables that can be easily measured from space. Calculation of ESI is typically dependent on complex parameterizations of evapotranspiration and knowledge of the surface energy budget, with limited utility if all input variables are not known or if measurement error is large. Evaluations of our proposed ESI parameterization using both a planetary boundary layer model and eddy covariance tower data from FLUXNET 2015 indicate that the model is a good approximation of ESI, particularly when the power law parameter is calculated by ecotype. Our analysis suggests that biosphere- atmosphere interactions, which directly impact vapor pressure deficit and relative humidity, can provide a simple but effective measure of soil moisture deficits.

# Notebooks 
Descriptions of the various Jupyter Notebooks used in the project

## data_wrangling.ipynb 
Using hourly [FLUXNET](https://fluxnet.org/) data, compute monthly means and add IGBP information from the Biological, Ancillary, Disturbance, and Metadata (BADM) file as a column to the table. Output a csv file for use in other notebooks.
#### *Input*
 - FLUXNET 2015 SUBSET hourly data for all tower sites
 - BADM file from FLUXNET 2015
#### *Output*
 - csv file of time averaged data
#### *Notebook dependencies*
 - pandas
 - numpy
 - glob 
 - sys 
 - os
 
## utils.ipynb
With FLUXNET 2015 monthly averaged data, the evaporative stress index (ESI) is computed using the Preistly Taylor equation. The exponential relationship between ESI and relative humidity is then used in conjunction with a simple parameterization of the stress index (detailed below and in [Fisher et al. 2008](https://www.sciencedirect.com/science/article/abs/pii/S0034425707003938)) in order to compute a scalar Beta value describing the curve. Functions defined in this notebook are used in other notebooks to solve for β and evaluate the model.
#### *Input*
 - FLUXNET 2015 monthly csv from data_wrangling.ipynb
#### *Output*
 - Figures
#### *Notebook dependencies*
 - pandas
 - numpy
 - os
 - scipy 
 - matplotlib.pyplot

## FLUXNET_ecotype_analysis.ipynb
Evaluate the evaporative stress index parameterization at different IGBP ecotypes using FLUXNET 2015 data. β values are calculated for each ecotype, and data is grouped by ecotype depending on the success of the parameterization in describing the Priestly Taylor evaporative stress index.
#### *Input*
 - FLUXNET 2015 monthly csv from data_wrangling.ipynb
#### *Output*
 - Figures
#### Notebook dependencies
 - pandas
 - numpy
 - os
 - scipy 
 - matplotlib.pyplot
 - sklearn
 - seaborn
 - nbimporter
 
## PBL_model_analysis.ipynb
This notebook analyzes output from Prof Dennis Baldocchi's MATLAB code for a growing and humidifying planetary boundary layer. See thesis for more information on the model.  
#### *Input*
 - None
#### *Output*
 - Figures
#### *Notebook dependencies*
 - pandas
 - numpy
 - matplotlib.pyplot
 - sklearn
 - seaborn
