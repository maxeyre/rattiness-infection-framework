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
This contains a set of parameters for controlling the rattiness model you wish to fit, each column is used as follows:
- `rat`: rattiness covariates you wish to include.
- `rat_par0`: starting values for rattiness parameters (regression coefficients, phi). Set to `FALSE` to fit model with starting values equal to zero.
- `with_Ui`: set to `TRUE` or `FALSE` to include/not include location-level iid random effects.
- `psi0`: starting value for $\psi$.

File: [control_R-rat_explore_R_H_spat_withnugg.csv](https://github.com/maxeyre/rattiness-infection-framework/blob/main/Data/control_R-rat_explore_R_H_spat_withnugg.csv)

#### 5. Model control parameters for joint rattiness-infection model
This contains a set of parameters for controlling the joint rattiness-infection model you wish to fit, each column is used as follows:
- `human`: human covariates you wish to include.
- `rat`: rattiness covariates you wish to include.
- `rat_par0`: starting values for rattiness parameters (regression coefficients, phi). Set to `FALSE` to fit model with starting values equal to zero.
- `par0.human`: starting values for human regression coefficients. Set to `TRUE` to start with values estimated by `glm()` or `FALSE` to fit model with starting values equal to zero.
- `multi.xi.on`: set to `TRUE` to fit a model with multiple $\xi$ parameters to test for interactions or `FALSE` to fit a model with a single $\xi$ parameter.
- `xi.var`: the variable (factor) you wish to use to define different $\xi$ parameter levels.
- `xi0`: starting values for xi (include as many as `xi.var` levels).
- `with_Ui`: set to `TRUE` or `FALSE` to include/not include normally-distributed location-level iid random effects in rattiness.
- `psi0`: starting value for the scale of spatial correlation of rattiness, $\psi$.
- `with_human_S`: set to `TRUE` or `FALSE` to include/not include a spatial Gaussian process in the human data. *(This was not included in the analysis in this paper because there was no residual spatial correlation in the human data)*
- `omega2.0`: starting value for the variance of the spatial Gaussian process in the human data. *(This was not included in the analysis in this paper)*
- `zeta0`: starting value for the scale of spatial correlation of spatial Gaussian process in the human data *(This was not included in the analysis in this paper)*
- `with_human_N`: set to `TRUE` or `FALSE` to include/not include normally-distributed location-level iid random effects in the human data.
- `omega2_nugg0`: starting value for $\omega^2$, the variance of the iid random effects in human data.

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
A description of all of the functions used in this analysis are given here to enable other users to use the rattiness-infection framework for their own data. If you would like to adapt the 

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
