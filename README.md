# Linking rattiness, geography and environmental degradation to spillover Leptospira infections in marginalised urban settings

## Article
The manuscript will be made available when published.

## Data
#### 1. Rat data
File: [1-PdLRattinessData.csv](https://github.com/maxeyre/Rattiness-1/blob/master/Data/1-PdLRattinessData.csv) 

The cleaned dataset for the cross-sectional survey for rat signs, traps and track plates in Pau da Lima community, Salvador, Brazil, which was collected between October and December 2014. It is georeferenced and includes all relevant covariate values used in the model described in the article. A codebook for variables can be found [here](https://github.com/maxeyre/Rattiness-1/blob/master/codebook.md). This dataset was collected by and belongs to the [Instituto de Saúde Coletiva, Universidade Federal da Bahia (ISC-UFBA)](http://www.isc.ufba.br/). If you use this dataset please could you cite ISC-UFBA and this published article (see license at bottom of page).

#### 2. Human data
File: [1-PdLHumanData.csv](https://github.com/maxeyre/Rattiness-1/blob/master/Data/1-PdLHumanData.csv) 

The cleaned dataset for the human cohort study in Pau da Lima community, Salvador, Brazil. It is georeferenced and includes all relevant covariate values used in the model described in the article. A codebook for variables can be found [here](https://github.com/maxeyre/Rattiness-1/blob/master/codebook.md). This dataset was collected by and belongs to the [Instituto de Saúde Coletiva, Universidade Federal da Bahia (ISC-UFBA)](http://www.isc.ufba.br/). If you use this dataset please could you cite ISC-UFBA and this published article (see license at bottom of page).

#### 3. Prediction grid
File: [2-PdLPredictionGrid.csv](https://github.com/maxeyre/Rattiness-1/blob/master/Data/2-PdLPredictionGrid.csv)

A 2.5m by 2.5m prediction grid with covariate values within the Pau da Lima study area.

## Software
The statistical software developed for this project and used for the analysis is divided into four main scripts and six additional scripts for fitting two-indices models and making predictions. Each script includes a description of its function within the R script. All scripts are contained [here](https://github.com/maxeyre/Rattiness-1/tree/master/Scripts) and are as follows:

#### 1. Exploratory analysis - Explore covariate relationships
File: [1-ExploreCovariates.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/1-ExploreCovariates.R)

#### 2. Exploratory analysis - Fit non-spatial models with covariates, LRT tests for each index, test for residual spatial correlation
File: [2-LRTTestSpatialCorr.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/2-LRTTestSpatialCorr.R)

#### 3. Full model - Fit multivariate geostatistical model
File: [3-FitFullModel.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/3-FitFullModel.R)

#### 4. Prediction - Make *rattiness* predictions at unobserved locations
File: [4-Prediction.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/4-Prediction.R)

#### 5. Two-indices models - Fitting and prediction
Folder: [Two-indices-models](https://github.com/maxeyre/Rattiness-1/tree/master/Scripts/Two-indices-models)

## Contact
If you have any questions about the project, data or software please contact max.eyre@lstmed.ac.uk.

## Licenses
#### Data:
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

#### Software: 
[License](https://github.com/maxeyre/Rattiness-1/blob/master/LICENSE)
