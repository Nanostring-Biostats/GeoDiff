#### Specs for fitPoisBG:
When groupvar is not provided,   
1. The function outputs a GeoMx S4 class with length same as length of ROIs, sizefact, in phenoData.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L20
2. The function outputs a GeoMx S4 class with length same as length of negative probes, featfact, in featureData. The value is NA for non-negative probes.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L23
3. If size_scale="first", sizefact[1]=1  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L30
4. If size_scale="sum", sum(sizefact)=1  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L34

When groupvar is provided but not found in the phenodata,   
1. The function returns an error saying this groupvar is not found in the S4 object.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L49 

When groupvar is provided and found in the phenodata, but the it only has one unique value,   
1. The function returns a warning message saying that the groupvar has only one value  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L63

When groupvar is provided and found in the phenodata with more than one unique value,     
1. The function outputs a GeoMx S4 class with length same as length of ROIs, sizefact_sp, in phenoData.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L95
2. The function outputs a GeoMx S4 class with length same as length of negative probes, featfact_sp, in featureData for each unique slide value. The value is NA for non-negative probes.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L91
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L99
3. If size_scale="first", sizefact[1]=1  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L115
4. If size_scale="sum", sum(sizefact)=1  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisBG.R#L119


#### Specs for diagPoisBG:
1. When generate_ppplot is TRUE, a figure is generated, when FALSE, no figure is generated  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-diagPoisBG.R#L19
2. when padj=FALSE, each element of sum of lowtail_prob and uptail_prob in the assay slot named lowtail_prob and uptail_prob equals to 1  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-diagPoisBG.R#L22
3. It returns an error without running fitPoisBG.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-diagPoisBG.R#L87
4. It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-diagPoisBG.R#L101


#### Specs for QuanRange:
1. The function outputs a GeoMx S4 class with length same as length of sample IDs (rownames) in phenoData for each probs input. The colname is the input prob.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-QuanRange.R#L43
2. It returns an error without running fitPoisBG.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-QuanRange.R#L78
3. It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-QuanRange.R#L91


#### Specs for BGScoreTest:
When split is FALSE,  
1. The function outputs a GeoMx S4 class with p values in featureData with length same as length of probes. The p value is NA for negative targets.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L41
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L87
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L132
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L177
2. The function outputs a GeoMx S4 class with score values in featureData with length same as length of targets. The score value is NA for negative probes.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L52
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L97
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L142
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L187
3. All p values are between 0 and 1 (inclusive) for non-negative features.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L62
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L107
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L152
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L197
4. The length of non-NA p values is equal to the number of non-negative probes.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L67
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L112
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L157
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L202
5. The length of non-NA scores values is equal to the number of non-negative probes.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L70
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L115
6. The order of pvalues is the same as scores.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L73
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L118
7. It returns an error without running fitPoisBG.    
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L508
8. It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L534

When split is TRUE,  
1. The function outputs a GeoMx S4 class with p values in featureData with length same as length of targets for each unique slide value. The p value is NA for negative probes.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L246
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L309
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L373
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L437
2. The function outputs a GeoMx S4 class with score values in featureData with length same as length of targets for each unique slide value. The score value is NA for negative probes.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L260
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L323
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L387
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L450
3. All p values are between 0 and 1 (inclusive) for non-negative features.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L272
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L335
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L399
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L462
4. The length of non-NA p values is equal to the number of non-negative probes.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L278
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L341
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L405
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L468
5. The length of non-NA scores values is equal to the number of non-negative probes.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L291
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L354
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L418
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L481
6. The order of pvalues is the same as scores.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L297
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L360
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L424
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-scoretest.R#L487


#### Specs for fitNBth:
1. Without providing values for features_high, sizefact_BG, threshold_start, the function returns the same value   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L68
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L171
2. The function outputs a GeoMx S4 class with para0 in the experimentData as "NA".   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L74
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L178
3. The function outputs a GeoMx S4 class with para, a matrix of estimated parameters, in the featureData. This matrix has feature_high_fitNBth in columns(same as features_high) and parameters(signal, r) in columns.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L80
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L184
4. The function outputs sizefact_fitNBth in the phenoData, which is positive, same length as sizefact_BG   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L92
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L196
5. The function outputs threshold in the experimentData. When threshold_fix=TRUE, threshold in the output is the same as threshold_start.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L98
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBth.R#L202


#### Specs for fitPoisthNorm:
When split is FALSE,  
1. Without providing values for ROIs_high, features_high, features_all, sizefact_start, sizefact_BG, the function returns the same value   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L128
2. User need to set iterations =2 for now   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L77
3. The function outputs a GeoMx S4 class with para0, matrix of estimated parameters, in the featureData.   
  This matrix para0 has the following structure:   
  3.1) 1 row for each feature (row name). If a feature is not in features_high, all columns will be NA.   
  3.2) n+1 columns labeled var1, var2, ..., var<n>, var<n+1> where n is the length of ROIs_high elements.   
  3.3) the n columns will have log2 expression (if feature is in features_high) or NA (otherwise).  
  3.4) the n+1th column contains the threshold for each feature in features_high and NA otherwise.    
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L134
4. The function outputs a GeoMx S4 class with para, matrix of estimated parameters, in the featureData. This matrix para has the following structure:  
  4.1) 1 row for each feature (row name) which is equal to the length of features_all  
  4.2) n+1 columns labeled var1, var2, ..., var<n>, var<n+1> where n is the length of ROIs_high elements.  
  4.3) the n columns will have log2 expression.  
  4.4) the n+1th column contains the threshold for each feature in features_all.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L153
5. The function outputs a column called conv0 in featureData, with values in [NA, 0] and length of 0s are the same as the length of features_high.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L169
6. The function outputs a column called conv in featureData, length same as features_all, and has values [NA, 0, 1]. The length of NAs equals the number of negative probes.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L178
7. when confac=0 and prior_type="contrast", preci1_norm value in  
   experimentData will be single value repeated over an n-by-n matrix (where n is the length of ROI_high).  
   This single value is equivalent, within 10 digits, to preci1con/n^2.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L189
8. It returns an error without running fitPoisBG or fitNBth.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L222
9. It returns an error if split is TRUE but no corresponding fitPoisBG is called previously.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L267

When split is TRUE,  
1. Without providing values for ROIs_high, features_high, features_all, sizefact_start, sizefact_BG, the function returns the same value  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L385
2. Given a GeoMx S4 object, fitPoisthNorm_sp runs the Poisson model-based normalization and log2 transformation on each element in "groupvar" individually. As such, the results for a given grouping/facet of the data should match the fitPoisthNorm results when an object is subset down to a single slide. Specifically, the following should be true:  
test:  
  2.1) For a given element in groupvar, the corresponding column in the "threshold0" matrix, which is within featureData, should be identical to the single-patient case's fetureData's para0[,n+1]th column.  
  2.2) For a given element in groupvar, the corresponding column in the "threshold" matrix, which is within featureData, should be identical to the single-patient case's fetureData's para[,n+1]th column.  
  2.3) For a given element in groupvar, the normalized matrix "normmat0_sp", in the assayData slot, should be identical to that element's "normmat0" matrix, also in the assayData slot, for all samples within that element. In other words, the matrix within the "single slide" results (normmat0) should be a subset of the "multiple slide" results (normmat_sp).  
  2.4) For a given element in groupvar, the normaized matrix "normmat_sp", in the assayData slot, should be identical to that element's "normmat" matrix, also in the assayData slot, for all samples within that element. In other words, the matrix within the "single slide" results (normmat) should be a subset of the "multiple slide" results (normmat_sp).  
  2.5) For a given element in groupvar, the vector of sizefact, located in phenoData's sizefact_norm column, is identical to that element's sizefact_norm vector from running fitPoisthNorm (i.e., single grouping case).  
  2.6) For a given element in groupvar, the vector of sizefact0, located in phenoData's sizefact_norm column, is identical to that element's sizefact_norm0 vector from running fitPoisthNorm (i.e., single grouping case).  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitPoisthNorm.R#L413

  
#### Specs for aggreprobe:
1. The function shall aggregate the probes depending on the argument provided.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L21
2. The function returns a GeoMxSet object when given a GeoMxSet object as input.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L33
3. For negative probes, the expression matrix and the target feature data available prior to collapsing shall match after collapsing except TargetName.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L38
4. TargetName for negative probes shall be updated to probe IDs after collapsing.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L48
5. For the non-negative probes, the subset of probes selected by the use method will be aggregated by sum into one target count.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L53
6. The resulting object shall have the same size of feature names as negative probe names plus non-negative target names.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L65
7. Single probe targets shall be returned without aggregation
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-aggreprobe.R#L73

#### Specs for fitNBthDE:
1.The function outputs para0, a matrix of estimated parameters in iter=1. This matrix has features_high in the columns and parameters(regression coefficients, threshold, r) in the rows.  Both threshold and r are positive.    
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthDE.R#L75
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthDE.R#L238
2.The function outputs para, a matrix of estimated parameters in iter=2. This matrix has features_all in the columns and parameters(regression coefficients, threshold, r) in the rows.  Both threshold and r are positive.   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthDE.R#L114
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthDE.R#L276
3.The function outputs sizefact, a vector of size factors, when sizefactrec=FALSE, sizefact is the same as sizefact_start.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthDE.R#L128
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthDE.R#L290


#### Specs for fitNBthmDE:
1. The function outputs para, a matrix of estimated parameters. This matrix has features_all in the columns and parameters(regression coefficients, threshold, r) in the rows. Both threshold and r are positive.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthmDE.R#L89
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-fitNBthmDE.R#L216


#### Specs for coefNBth:
1. when fullpara=TRUE, the output parameters should be regression coefficients, threshold and r in a list. Both threshold and r are positive.  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-coefNBth.R#L22
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-coefNBth.R#L57
2. when fullpara=FALSE, the output parameters should be regression coefficients only in a list  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-coefNBth.R#L30
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-coefNBth.R#L65


#### Specs for contrastNBth:
1. The function takes in a DE model as an input from fitNBthDE or fitNBthmDE  
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-contrastNBth.R#L19
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-contrastNBth.R#L59
2. The user input test:statistical test, choose from c("two-sided", ">", "<")   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-contrastNBth.R#L20
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-contrastNBth.R#L60
3. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-contrastNBth.R#L30
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-contrastNBth.R#L70


#### Specs for DENBth:
1. In the output data.frame, the DE table should not have NA if NAto1=TRUE   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-DENBth.R#L22
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-DENBth.R#L54
2. In the output data.frame, the DE table has a column padj when adjp=TRUE   
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-DENBth.R#L25
test: https://github.com/Nanostring-Biostats/GeoDiff/blob/eef13efc1636fd86e3b833c26967d3fb2c350396/tests/testthat/test-DENBth.R#L57