

#### Reqs for fitPoisBG:
- The user input object, a GeoMx S4 class
- The user input the group variable name for slide id of ROIs: groupvar
- The user input the size_scale from "sum", "first" for sizefact (default "sum")
- The user input iteration tolerance: tol (default 1e-3)
- The iteration number: iterations (default 10)
- The function outputs a GeoMx S4 class 
  - The function outputs sample size factors vector in phenoData: sizefact or size_fact_sp if if multiple slides exist
  - The function outputs feature factor vector in featureData: featfact or featfeact_sp_XXX if multiple slides exist (XXX is the slide name)
  - the group variable name for slide id of ROIs in experimentData: fitPoisBG_sp_var if multiple slides exist
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-fitpoisbg


#### Reqs for diagPoisBG:
- The user input object, a GeoMx S4 class
- The user input an indicator variable on whether using results from multiple slides: split
- The user input whether to adjust p value: padj (default TRUE)
- The user input adjust p value method: padj_method (default "BH") 
- The user input p value: cutoff (default 1e-6)
- The user input whether to generate ppplot: generate_ppplot(default TRUE) 
- The function generates ppplot if generate_ppplot=TRUE
- The function outputs a GeoMx S4 class
  - matrix of lower tail probability in assayData slot called lowtail_prob : lowtail_prob
  - matrix of upper tail probability in assayData slot called uptail_prob: uptail_prob 
  - the dispersion parameter in experimentData named disper or disper_sp if multiple is TRUE: disper or disper_sp
  - matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called low_outlier : low_outlier
  - matrix of outlier indicator (Yes: 1; No: 0) in assayData slot called up_outlier : up_outlier
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-diagpoisbg


#### Reqs for QuanRange
- The user input a GeoMx S4 class
- The user input an indicator variable on whether using results from multiple slides: split
- The user input vector of probabilities: probs (default FALSE)
- The function outputs a GeoMx S4 class with column names same as probs:quanrange in phenoData
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-quanrange


#### Reqs for BGScoreTest:
- The user input a GeoMx S4 class
- The user input an indicator variable on whether using results from multiple slides: split
- The user input adjustment factor: adj (default 1)
- The user input whether to use outlier detection and remove outlier: removeoutlier (default FALSE)
- The user input whether to implement the score test with the assumption the feature factor follows gamma distribution: useprior (default FALSE)
- The function outputs a GeoMx S4 class if split is FALSE. 
  - vector of p values in featureData:pvalues
  - vector of Background score test statistics in featureData: scores
- The function outputs a GeoMx S4 class if split is TRUE. 
  - matrix of p values with column dimension equals to slide number in featureData: pvalues_XXX
  - matrix of Background score test statistics with column dimension equals to slide number in featureData: scores_XXX
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-bgscoretest


#### Reqs for fitNBth
- The user input a GeoMx S4 class
- The user input vector of high features names: features_high (If NULL, the default value is calculated within the function)
- The user input size factor for background: sizefact_BG (If NULL, the default value is calculated within the function)
- The user input starting size factor: sizefact_start(default=sizefactor_BG)
- The user input the size_scale from "sum", "first" for sizefact (default "sum")
- The user input starting threshold: threshold_start (If NULL, the default value is calculated within the function)
- The user input whether to fix threshold: threshold_fix(default=FALSE)
- The user input iteration tolerance: tol (default=1e-3)
- The user input iteration number: iterations (default=5)
- The user input start_para: starting values for parameter estimation (default=c(threshold_start, 1))
- The user input lower_sizefact: lower limit for sizefact (default=0)
- The user input the lower end of threshold: lower_threshold(default=threshold_start/5)
- The function outputs a GeoMx S4 class
  - para0 in experimentData called para0 (=NA)
  - para in featureData called para, matrix of estimated parameters, features in rows and parameters(signal, r) in columns.
  - sizefact in phenoData called sizefact_fitNBth, estimated sizefact
  - preci1 in experimentData call preci1 (=NA)
  - conv0 in experimentData call conv0 (=NA)
  - conv in experimentData call conv(=NA)
  - Im in experimentData call Im(=NA)
  - features_high in featureData called feature_high_fitNBth, same as the input features_high
  - features_all in experimentData called features_all(=NA)
  - threshold in experimentData called threshold
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-fitnbth


#### Reqs for fitPoisthNorm  
- The user input a GeoMx S4 class
- The user input an indicator variable on whether using results from multiple slides: split
- The user input vector of high ROIs names: ROIs_high (If NULL, the default value is calculated within the function)
- The user input vector of high features names: features_high (If NULL, the default value is calculated within the function)
- The user input vector of all features names need to be fitted: features_all (If NULL, the default value is calculated within the function)
- The user input size factor for background: sizefact_BG (If NULL, the default value is calculated within the function)
- The user input starting size factor: sizefact_start (If NULL, the default value is calculated within the function)
- The user input threshold_mean: threshold_mean  
- The user input precision matrix for threshold: preci2  
- The user input iteration number: iterations=1 or 2, default=2
- The user input prior type for preci1: prior_type from "equal" or "contrast", default="contrast"
- The user input XXX: sizefactrec, default = TRUE
- The user input the size_scale from "sum", "first" for sizefact (default "sum")
- The user input XXX: sizescalebythreshold, default = TRUE
- The user input robust covariance: covrob, default=FALSE
- The user input constant term in specifying precision matrix 1: preci1con(default=1/25) 
- The user input cutoff term in calculating precision matrix 1 (default=15) 
- The user input factor for contrast in precision matrix 1:confac (default=1) 
- The user input whether to calculate hessian: calhes (default=FALSE)
- When split is FALSE, the function outputs a GeoMx S4 class
  - para0 in featureData called para0, matrix of estimated parameters for iter=1, features in columns and parameters(log2 expression, threshold) in rows.
  - para in featureData called para, matrix of estimated parameters for iter=2, features in columns and parameters(log2 expression, threshold) in rows.
  - normmat0 in assayData slot called normmat0, matrix of log2 expression for iter=1, features in columns and log2 expression in rows.
  - normmat in assayData slot called normmat, matrix of log2 expression for iter=2, features in columns and log2 expression in rows. 
  - sizefact in phenoData called sizefact_norm, estimated sizefact
  - sizefact0 in phenoData called sizefact0_norm, estimated sizefact in iter=1
  - preci1 in featureData called preci1_norm, precision matrix 1
  - conv0 in featureData called conv0, vector of convergence for iter=1, 0 converged, 1 not converged 
  - conv in featureData called conv, vector of convergence for iter=2, 0 converged, 1 not converged
  - features_high in featureData called features_high , same as the input features_high
  - features_all in featureData called features_all, same as the input features_all
- When split is TRUE, the function outputs a GeoMx S4 class
  - threshold0 in featureData called threshold0, matrix of estimated threshold for iter=1, features in columns and threshold for different slides in rows.
  - threshold in featureData called threshold, matrix of estimated threshold for iter=2, features in columns and threshold for different slides in rows.
  - normmat0 in assayData slot called normmat0_sp, matrix of log2 expression for iter=1, features in columns and log2 expression in rows.
  - normmat in assayData slot called normmat_sp, matrix of log2 expression for iter=2, features in columns and log2 expression in rows. 
  - sizefact in phenoData called sizefact_norm_sp, estimated sizefact
  - sizefact0 in phenoData called sizefact0_norm_sp, estimated sizefact in iter=1
  - preci1 in experimentData called preci1_norm_sp, precision matrix 1
  - conv0 in featureData called conv0_sp_XX, vector of convergence for iter=1, 0 converged, 1 not converged, NA if features are not present. 
  - conv in featureData called conv_sp_XX, vector of convergence for iter=2, 0 converged, 1 not converged, NA if features are not present. 
  - features_high in featureData called features_high_sp, same as the input features_high
  - features_all in featureData called features_all_sp, same as the input features_all
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-fitpoisthnorm


#### Reqs for aggreprobe
- The user input a GeoMx S4 class
- The user input an indicator variable on whether using results from multiple slides: split
- The user input the method to determine outliers from score, cor, and both
- The function outputs a GeoMx S4 class 
  - the collapsed expression matrix by target in the expression slot 
  - the collapsed feature data by target in the featureData slot
  - probenum in featureData, numerical vector of probe numbers of targets
  - proberemained in featureData, the list of remaining probes of targets
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-aggreprobe


#### Reqs for fitNBthDE
- The user input a GeoMx S4 class 
- The user input form: model formula
- The user input an indicator variable on whether using results from multiple slides: split
- The user input features_high: features with high abudance to be run in iter=1
- The user input features_all: all features to be run in iter=2 
- The user input sizefact_start: initial value for size factors
- The user input sizefact_BG: size factor for background
- The user input threshold_mean: average threshold level
- The user input preci2: precision for the background
- The user input lower_threshold: lower limit for the threshold
- The user input prior_type: empirical bayes prior type, choose from c("equal", "contrast")
- The user input sizefactrec, whether to recalculate sizefact, default=TRUE
- The user input size_scale: how to scale size factor if sizefactrec=TRUE
- The user input sizescalebythreshold: whether to scale the size factor by the threshold_mean in the modeling, default=TRUE
- The user input iterations: how many iterations need to run to get final results, default=2,
                   the first iteration apply the model only on features_high and construct the prior then refit the model using                     this prior for all genes.
- The user input covrob: whether to use robust covariance in calculating covariance. default=FALSE
- The user input preci1con: constant for preci1
- The user input cutoff: cutoff for calculating the precision matrix for regression coefficients
- The user input confac: contrast factor in the precision matrix for regression coefficients
- The function outputs a list of following objects
  - design matrix: X = X
  - parameters estimated in iter 1: para0
  - parameters estimated in iter 2: para
  - size factor for signal: sizefact
  - size factor for background: sizefact0,
  - preci matrix for regression coefficients, preci1,
  - Information matrix: Im0,
  - Information matrix: Im,
  - vector of whether model has converged in iter=1: conv0, 0=converged, 1=not converged
  - vector whether model has converged in iter=2: conv, 0=converged, 1=not converged
  - features with high abundance to be run in iter=1: features_high
  - all features need to be run in iter=2: features_all
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-fitnbthde


#### Reqs for fitNBthmDE

- The user input a GeoMx S4 class
- The user input form: model formula
- The user input an indicator variable on whether using results from multiple slides: split
- The user input features_all: vector of features to be run.
- The user input sizefact size factor: size factor for signal 
- The user input sizefact_BG: size factor for background
- The user input preci1: precision matrix for regression coefficients
- The user input threshold_mean:  average background level
- The user input preci2 precision for the background
- The user input seed random number generator seed
- The user input sizescalebythreshold: whether to scale the size factor, default=TRUE
- The user input controlRandom: list of random effect control parameters
- The function outputs a list of following objects
  - X, design matrix for fixed effect
  - Z, design matrix for random effect
  - rt, random effect terms
  - para0, =NA
  - para, estimated parameters, including regression coefficients, r and threshold in rows and features in columns
  - sizefact, same as input sizefact
  - sizefact0, NA
  - preci1, input precision matrix for regression coefficients
  - Im0, NA
  - Im, Information matrix of parameters
  - conv0, NA
  - conv, vector of convergence, 0 converged, 1 not converged
  - features_high, NA
  - features_all, same as the input features_all
  - theta, list of estimated random effect parameters(for relative covariance matrix)
  - varcov, list of estimated variance covariance parameter estimation
  - MAP random effect
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-fitnbthmde


#### Reqs for coefNBth
- The user input object: DE model, output by fitNBthDE or fitNBthmDE
- The user input fullpara: fullpara whether to generate results on all parameters
- The function outputs a list of following objects
  - estimate, coefficients estimate
  - wald_stat, Wald test statistics
  - p_value, p value of Wald test
  - se, standard error
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-coefnbth


#### Reqs for contrastNBth
- The user input object: DE model, output by fitNBthDE or fitNBthmDE
- The user input test: statistical test, choose from c("two-sided", ">", "<")
- The user input method: contrasts methods, only matrix of contrast vector is allowed for now, default=diag(1,ncol(object$X)), i.e. testing the regression coefficients
- The user input baseline: testing baseline, default=0.
- The function outputs a list of following objects
  - estimate, contrasts estimate
  - wald_stat, Wald test statistics
  - p_value, p value of Wald test
  - se, standard error
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-contrastnbth


#### Reqs for DENBth

- The user input object: inference list from coefNBth or contrastNBth
- The user input variable: variable contrasts need to construct on
- The user input NAto1 whether to replace NA in pvalue by 1
- The user input padj whether to adjust p value
- The user input padj_method p value adjustment method, default="BH"
- The function outputs a data frame with following columns 
  - log2FC, fold change in log scale
  - pvalue, unadjusted p values
  - adjp, adjusted p values
specifications: https://github.com/Nanostring-Biostats/GeoDiff/blob/main/specs.md#specs-for-denbth
