# Linking rattiness, geography and environmental degradation to spillover *Leptospira* infections in marginalised urban settings

## Article
The manuscript will be made available when published.

## Data
A codebook for all data can be found [here](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/codebook.md).
#### 1. Rat data
File: [rat_data.csv](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/rat_data.csv) 

The cleaned dataset for the cross-sectional survey in Pau da Lima community, Salvador, Brazil. It contains information on the five rat abundance indices (traps, track plates, number of burrows, presence of faeces and presence of trails), and was collected between October and December 2014. It is georeferenced and includes all relevant covariate values used in the model described in the article. This dataset was collected by and belongs to the [Instituto de Saúde Coletiva, Universidade Federal da Bahia (ISC-UFBA)](http://www.isc.ufba.br/). If you use this dataset please could you cite ISC-UFBA and this published article (see license at bottom of page).

#### 2. Human data
File: [human_data.csv](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/human_data.csv) 

The cleaned dataset for the human cohort study in Pau da Lima community, Salvador, Brazil. It is georeferenced and includes all relevant covariate values used in the model described in the article. This dataset was collected by and belongs to the [Instituto de Saúde Coletiva, Universidade Federal da Bahia (ISC-UFBA)](http://www.isc.ufba.br/). If you use this dataset please could you cite ISC-UFBA and this published article (see license at bottom of page).

#### 3. Prediction grid
File: [prediction_grid.csv](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/prediction_grid.csv)

A 2.5m by 2.5m prediction grid with covariate values within the Pau da Lima study area.

#### 4. Model control parameters for rattiness model
File: [control_R-rat_explore_R_H_spat_withnugg.csv](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/control_R-rat_explore_R_H_spat_withnugg.csv)

#### 5. Model control parameters for joint rattiness-infection model
File: [control_R1_3Xi.csv](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/control_R1_3Xi.csv)

## Analysis scripts
R scripts to conduct the analysis described in this article are available [here](https://github.com/maxeyre/rattiness-infection-framework/tree/main/Scripts). Each script includes a brief description of its function.

#### 1. Exploratory analysis of rat data
File: [1-rat-explore.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/1-rat-explore.R)

#### 1. Exploratory analysis of human data
File: [2-human-explore.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/2-human-explore.R)

#### 3. Fit the joint rattiness-infection model and bootstrap for CIs
File: [3-fit-joint-model.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/3-fit-joint-model.R)

#### 4. Spatial prediction for rattiness and probability of infection
File: [4-spatial-prediction.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/4-spatial-prediction.R)

## Functions

#### 1. Exploratory analysis of rat data
File: [1-rat-explore.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/1-rat-explore.R)

#### 1. Exploratory analysis of human data
File: [2-human-explore.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/2-human-explore.R)

#### 3. Fit the joint rattiness-infection model and bootstrap for CIs
File: [3-fit-joint-model.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/3-fit-joint-model.R)

#### 4. Spatial prediction for rattiness and probability of infection
File: [4-spatial-prediction.R](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Scripts/4-spatial-prediction.R)


## Contact
If you have any questions about the project, data or software please contact max.eyre@lstmed.ac.uk.

## Licenses
#### Data:
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

#### Software: 
[License](https://github.com/maxeyre/rattiness-infection-framework/blob/main/LICENSE)
