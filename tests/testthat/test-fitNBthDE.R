test_that("fitNBthDE produces desired results, CTA", {
  
  #### Specs for fitNBth
  # 1 The function outputs para0, a matrix of estimated parameters in
  #   iter=1. This matrix has features_high in the columns and
  #   parameters(regression coefficients, threshold, r) in the rows.
  #   Both threshold and r are positive.
  # 2 The function outputs para, a matrix of estimated parameters in
  #   iter=2. This matrix has features_all in the columns and
  #   parameters(regression coefficients, threshold, r) in the rows.
  #   Both threshold and r are positive.
  # 3 The function outputs sizefact, a vector of size factors,
  #   when sizescalebythreshold=FALSE, sizefact is the same as
  #   sizefact_start.
  
  library(dplyr)
  ### Initializing CTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions
  data("demoData") # for tests of structure of demoData itself, see test-scoretest.R
  
  # susbet samples
  demoData <- demoData[, c(1:5, 33:37)]
  set.seed(413)
  demoData <- fitPoisBG(demoData, size_scale = "sum")
  demoData <- aggreprobe(demoData, use = "cor")
  demoData <- BGScoreTest(demoData)
  demoData$slidename <- substr(demoData[["slide name"]], 12, 17)
  thmean <- 1 * mean(fData(demoData)$featfact, na.rm = TRUE)
  demo_pos <- demoData[which(!fData(demoData)$CodeClass == "Negative"), ]
  demo_neg <- demoData[which(fData(demoData)$CodeClass == "Negative"), ]
  sc1_scores <- fData(demo_pos)[, "scores"]
  names(sc1_scores) <- fData(demo_pos)[, "TargetName"]
  features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) &
                      (sc1_scores < quantile(sc1_scores, probs = 0.95))) |>
    which() |>
    names()
  demoData <- fitNBth(demoData,
                      features_high = features_high,
                      sizefact_BG = demo_neg$sizefact,
                      threshold_start = thmean,
                      iterations = 5,
                      start_para = c(200, 1),
                      lower_sizefact = 0,
                      lower_threshold = 100,
                      tol = 1e-8)
  ROIs_high <- sampleNames(demoData)[which(demoData$sizefact_fitNBth * thmean > 2)]
  features_all <- rownames(demo_pos)
  
  pData(demoData)$group <- c(rep(1, 5), rep(2, 5))
  
  ### Case 1:
  
  features_high <- features_all
  
  NBthDEmod1 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = demoData,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = demoData[, ROIs_high][["sizefact_fitNBth"]],
    sizefact_BG = demoData[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    iterations = 1
  )
  
  # 1 The function outputs para0,...
  #   a matrix of estimated parameters in iter=1.
  expect_true(all(features_high == features_all))
  expect_true("para0" %in% names(NBthDEmod1))
  expect_true(inherits(NBthDEmod1[["para0"]], "matrix"))
  #   This matrix has features_high in the columns
  para0 <- NBthDEmod1[["para0"]]
  expect_true(ncol(para0) == length(features_high))
  expect_true(all(colnames(para0) == features_high))
  #   and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "group", "r", "threshold") == rownames(para0)))
  #   Both threshold and r are positive.
  expect_true(all(para0["r", ] > 0))
  expect_true(all(para0["threshold", ] > 0))
  
  ### Case 2:
  
  features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) &
                      (sc1_scores < quantile(sc1_scores, probs = 0.95))) |>
    which() |>
    names()
  
  NBthDEmod2 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = demoData,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = demoData[, ROIs_high][["sizefact_fitNBth"]],
    sizefact_BG = demoData[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    iterations = 2
  )
  
  # 2 The function outputs para,...
  #   a matrix of estimated parameters in iter=2.
  expect_true("para" %in% names(NBthDEmod2))
  expect_true(inherits(NBthDEmod2[["para"]], "matrix"))
  #   This matrix has features_all in the columns
  para <- NBthDEmod2[["para"]]
  expect_true(ncol(para) == length(features_all))
  expect_true(all(colnames(para) == features_all))
  #   and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "group", "r", "threshold") == rownames(para)))
  #   Both threshold and r are positive.
  expect_true(all(para["r", ] > 0))
  expect_true(all(para["threshold", ] > 0))
  
  # 3 The function outputs sizefact,...
  #   a vector of size factors,
  expect_true("sizefact" %in% names(NBthDEmod2))
  expect_true(is.vector(NBthDEmod2[["sizefact"]]))
  #   when sizescalebythreshold=FALSE,...
  
  ### Case 3:
  
  sizefact_start = demoData[, ROIs_high][["sizefact_fitNBth"]]
  set.seed(123)
  NBthDEmod3 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = demoData,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = sizefact_start,
    sizefact_BG = demoData[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    sizefactrec = FALSE
  )
  
  
  #   sizefact is the same as sizefact_start.
  expect_true("sizefact" %in% names(NBthDEmod3))
  expect_true(is.vector(NBthDEmod3[["sizefact"]]))
  expect_true(length(NBthDEmod3[["sizefact"]]) == length(sizefact_start))
  expect_true(all(NBthDEmod3[["sizefact"]] == sizefact_start))
  
  
})

test_that("fitNBthDE works with dgCMatrix format",{
  data("demoData")
  demoData <- demoData[, c(1:5, 33:37)]
  set.seed(413)
  demoData <- fitPoisBG(demoData, size_scale = "sum")
  demoData <- aggreprobe(demoData, use = "cor")
  negdat <- demoData[which(Biobase::fData(demoData)$CodeClass == "Negative"), ]
  countmat <- Biobase::exprs(negdat)
  countmat = as(countmat, "dgCMatrix")
  result <- fitPoisBG(
    object = countmat,
    iterations = 10,
    tol = 1e-3,
    size_scale = "sum")
  demoData[["sizefact"]] <- result$sizefact[Biobase::sampleNames(demoData)]
  Biobase::fData(demoData)[["featfact"]] <- NA
  Biobase::fData(demoData)[["featfact"]][match(names(result$featfact), Biobase::featureNames(demoData), nomatch = 0)] <- result$featfact
  # Case 1: adjustment factor 5, no outlier removal, no prior
  demoData <- BGScoreTest(demoData,
                          adj = 5,
                          removeoutlier = FALSE, useprior = FALSE)
  demoData$slidename <- substr(demoData[["slide name"]], 12, 17)
  thmean <- 1 * mean(fData(demoData)$featfact, na.rm = TRUE)
  demo_pos <- demoData[which(!fData(demoData)$CodeClass == "Negative"), ]
  demo_neg <- demoData[which(fData(demoData)$CodeClass == "Negative"), ]
  sc1_scores <- fData(demo_pos)[, "scores"]
  names(sc1_scores) <- fData(demo_pos)[, "TargetName"]
  features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) &
                      (sc1_scores < quantile(sc1_scores, probs = 0.95))) |>
    which() |>
    names()
  demoData <- fitNBth(demoData,
                      features_high = features_high,
                      sizefact_BG = demo_neg$sizefact,
                      threshold_start = thmean,
                      iterations = 5,
                      start_para = c(200, 1),
                      lower_sizefact = 0,
                      lower_threshold = 100,
                      tol = 1e-8)
  ROIs_high <- sampleNames(demoData)[which(demoData$sizefact_fitNBth * thmean > 2)]
  features_all <- rownames(demo_pos)
  
  pData(demoData)$group <- c(rep(1, 5), rep(2, 5))
  
  ### Case 1:
  
  features_high <- features_all
  
  NBthDEmod1 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = demoData,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = demoData[, ROIs_high][["sizefact_fitNBth"]],
    sizefact_BG = demoData[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    iterations = 1
  )
  
  # 1 The function outputs para0,...
  #   a matrix of estimated parameters in iter=1.
  expect_true(all(features_high == features_all))
  expect_true("para0" %in% names(NBthDEmod1))
  expect_true(inherits(NBthDEmod1[["para0"]], "matrix"))
  #   This matrix has features_high in the columns
  para0 <- NBthDEmod1[["para0"]]
  expect_true(ncol(para0) == length(features_high))
  expect_true(all(colnames(para0) == features_high))
  
})


test_that("fitNBthDE produces desired results, WTA", {
  
  ### Same overall workflow as above but with the WTA kidney dataset
  
  library(dplyr)
  ### Initializing WTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # Run data through (required) upstream functions
  data("kidney")
  set.seed(413)
  kidney <- kidney[, kidney$`slide name` %in% c("disease1B", "disease2B")]
  #kidney <- kidney[, c(1:10, 11:20)]
  kidney <- fitPoisBG(kidney, size_scale = "sum")
  kidney <- fitPoisBG(kidney, groupvar = "slide name", size_scale = "sum")
  all0probeidx <- which(rowSums(exprs(kidney))==0)
  kidney <- kidney[-all0probeidx, ]
  kidney <- aggreprobe(kidney, use = "cor")
  kidney <- BGScoreTest(kidney)
  # Negative and Non-Negative facets:
  kidney_neg <- kidney[which(fData(kidney)$CodeClass == "Negative"), ]
  kidney_pos <- kidney[-which(fData(kidney)$CodeClass == "Negative"), ]
  # feature factors per groupvar (i.e., slide name)
  featfact_sp <- fData(kidney_neg)[, grep("featfact_", fvarLabels(kidney_neg))]
  # scaling factors
  posdat <- Biobase::exprs(kidney_pos)
  gene_sum <- rowSums(posdat)
  features_high <- ((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95))) |>
    which() |>
    names()
  set.seed(123)
  genes_high <- sample(features_high, 1500) # subset
  thmean <- 1 * (featfact_sp |> colMeans())[1] # picks the first slide name's value
  
  kidney <- fitNBth(kidney,
                    features_high = genes_high,
                    sizefact_BG = kidney_neg$sizefact_sp,
                    threshold_start = thmean,
                    iterations = 5,
                    start_para = c(200, 1),
                    lower_sizefact = 0,
                    lower_threshold = 100,
                    #threshold_fix = FALSE, # default but calling it explicitly here
                    tol = 1e-8
  )
  
  ROIs_high <- sampleNames(kidney)[which(kidney$sizefact_fitNBth * thmean > 2)]
  features_all <- rownames(kidney_pos)
  
  pData(kidney)$group <- c(rep(1, 38), rep(2, 38))
  
  ### Case 1:
  
  features_high <- features_all
  
  NBthDEmod1 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = kidney,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = kidney[, ROIs_high][["sizefact_fitNBth"]],
    sizefact_BG = kidney[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    iterations = 1
  )
  
  # 1 The function outputs para0,...
  #   a matrix of estimated parameters in iter=1.
  expect_true(all(features_high == features_all))
  expect_true("para0" %in% names(NBthDEmod1))
  expect_true(inherits(NBthDEmod1[["para0"]], "matrix"))
  #   This matrix has features_high in the columns
  para0 <- NBthDEmod1[["para0"]]
  expect_true(ncol(para0) == length(features_high))
  expect_true(all(colnames(para0) == features_high))
  #   and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "group", "r", "threshold") == rownames(para0)))
  #   Both threshold and r are positive.
  expect_true(all(para0["r", ] > 0))
  expect_true(all(para0["threshold", ] > 0))
  
  ### Case 2:
  
  features_high <- ((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95))) |>
    which() |>
    names()
  
  NBthDEmod2 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = kidney,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = kidney[, ROIs_high][["sizefact_fitNBth"]],
    sizefact_BG = kidney[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    iterations = 2
  )
  
  # 2 The function outputs para,...
  #   a matrix of estimated parameters in iter=2.
  expect_true("para" %in% names(NBthDEmod2))
  expect_true(inherits(NBthDEmod2[["para"]], "matrix"))
  #   This matrix has features_all in the columns
  para <- NBthDEmod2[["para"]]
  expect_true(ncol(para) == length(features_all))
  expect_true(all(colnames(para) == features_all))
  #   and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "group", "r", "threshold") == rownames(para)))
  #   Both threshold and r are positive.
  expect_true(all(para["r", ] > 0))
  expect_true(all(para["threshold", ] > 0))
  
  # 3 The function outputs sizefact,...
  #   a vector of size factors,
  expect_true("sizefact" %in% names(NBthDEmod2))
  expect_true(is.vector(NBthDEmod1[["sizefact"]]))
  #   when sizescalebythreshold=FALSE,...
  
  ### Case 3:
  
  sizefact_start = kidney[, ROIs_high][["sizefact_fitNBth"]]
  set.seed(123)
  NBthDEmod3 <- fitNBthDE(
    form = ~group,
    split = FALSE,
    object = kidney,
    ROIs_high = ROIs_high,
    features_high = features_high,
    features_all = features_all,
    sizefact_start = kidney[, ROIs_high][["sizefact_fitNBth"]],
    sizefact_BG = kidney[, ROIs_high][["sizefact"]],
    preci2 = 10000,
    prior_type = "contrast",
    covrob = FALSE,
    preci1con = 1 / 25,
    sizescalebythreshold = TRUE,
    sizefactrec = FALSE
  )
  
  
  #   sizefact is the same as sizefact_start.
  expect_true("sizefact" %in% names(NBthDEmod3))
  expect_true(is.vector(NBthDEmod3[["sizefact"]]))
  expect_true(length(NBthDEmod3[["sizefact"]]) == length(sizefact_start))
  expect_true(all(NBthDEmod3[["sizefact"]] == sizefact_start))

  
  
})

test_that("coefNBth produces desired results from output of fitNBthmDE", {
 library(dplyr)
 ### Initializing CTA objects before running tests
 # Create temporary directory that will get destroyed after this block is executed.
 tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
 withr::local_dir(tmp_dir)
 # Run data through (required) upstream functions
 data("kidney")
 set.seed(413)
 kidney <- kidney[, kidney$`slide name` %in% c("disease1B", "disease2B")]
 #kidney <- kidney[, c(1:10, 11:20)]
 kidney <- fitPoisBG(kidney, size_scale = "sum")
 kidney <- fitPoisBG(kidney, groupvar = "slide name", size_scale = "sum")
 all0probeidx <- which(rowSums(exprs(kidney))==0)
 kidney <- kidney[-all0probeidx, ]
 kidney <- aggreprobe(kidney, use = "cor")
 kidney <- BGScoreTest(kidney)
 # Negative and Non-Negative facets:
 kidney_neg <- kidney[which(fData(kidney)$CodeClass == "Negative"), ]
 kidney_pos <- kidney[-which(fData(kidney)$CodeClass == "Negative"), ]
 # feature factors per groupvar (i.e., slide name)
 featfact_sp <- fData(kidney_neg)[, grep("featfact_", fvarLabels(kidney_neg))]
 # scaling factors
 posdat <- Biobase::exprs(kidney_pos)
 gene_sum <- rowSums(posdat)
 features_high <- ((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95))) |>
  which() |>
  names()
 set.seed(123)
 genes_high <- sample(features_high, 1500) # subset
 thmean <- 1 * (featfact_sp |> colMeans())[1] # picks the first slide name's value
 
 kidney <- fitNBth(kidney,
                   features_high = genes_high,
                   sizefact_BG = kidney_neg$sizefact_sp,
                   threshold_start = thmean,
                   iterations = 5,
                   start_para = c(200, 1),
                   lower_sizefact = 0,
                   lower_threshold = 100,
                   #threshold_fix = FALSE, # default but calling it explicitly here
                   tol = 1e-8
 )
 
 ROIs_high <- sampleNames(kidney)[which(kidney$sizefact_fitNBth * thmean > 2)]
 features_all <- rownames(kidney_pos)
 
 pData(kidney)$group <- c(rep(1, 38), rep(2, 38))
 
 ### Case 1:
 
 features_high <- features_all
 features_high <- ((gene_sum > quantile(gene_sum, probs = 0.5)) & (gene_sum < quantile(gene_sum, probs = 0.95))) |>
  which() |>
  names()
 
 NBthDEmod2 <- fitNBthDE(
  form = ~group,
  split = FALSE,
  object = kidney,
  ROIs_high = ROIs_high,
  features_high = features_high,
  features_all = features_all,
  sizefact_start = kidney[, ROIs_high][["sizefact_fitNBth"]],
  sizefact_BG = kidney[, ROIs_high][["sizefact"]],
  preci2 = 10000,
  prior_type = "contrast",
  covrob = FALSE,
  preci1con = 1 / 25,
  sizescalebythreshold = TRUE,
  iterations = 2
 )

 coefmfull <- coefNBth(NBthDEmod2, fullpara=TRUE)
 coefmreg <- coefNBth(NBthDEmod2, fullpara=FALSE)
 
 ## 1. when fullpara=TRUE, the output parameters should be regression coefficients, threshold and r in a list.
 ##    Both threshold and r are positive.
 
 
 expect_true(all(c(colnames(NBthDEmod2$X), c("r", "threshold")) == rownames(coefmfull$estimate)))
 expect_true(all(coefmfull$estimate["r", ] > 0))
 expect_true(all(coefmfull$estimate["threshold", ] > 0))
 
 ## 2. when fullpara=FALSE, the output parameters should be regression coefficients only in a list
 
 expect_true(all(colnames(NBthDEmod2$X) == rownames(coefmreg$estimate)))
 
 ## 3. The user input test:statistical test, choose from c("two-sided", ">", "<")
 ## 4. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1
 coeffull <- contrastNBth(NBthDEmod2)
 coeftest <- contrastNBth(NBthDEmod2, method=matrix(c(0,1), 2, 1), baseline=0)
 
 coeftestupper <- contrastNBth(NBthDEmod2, method=matrix(c(0,1), 2, 1), baseline=0, test = ">")
 
 coeftestlower <- contrastNBth(NBthDEmod2, method=matrix(c(0,1), 2, 1), baseline=0, test = "<")
 
 
 
 expect_true(all(coeffull$estimate[2,] == coeftest$estimate[1,]))
 
 ## 5. In the output list, the p values of '>' and '<' for the same variable/feature should add up to 1
 
 expect_equal(unname(coeftestlower$p_value+coeftestupper$p_value),
              rep(1, length(coeftestlower$p_value)))
})


test_that("coefNBth produces desired results from output of fitNBthmDE and parallel works properly", {
 library(dplyr)
 ### Initializing CTA objects before running tests
 # Create temporary directory that will get destroyed after this block is executed.
 tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
 withr::local_dir(tmp_dir)
 # Run data through (required) upstream functions
 library(Rfast)
 data("test_data")
 
 ## load data
 #library(GeoDiff)
 library(Matrix)
 library(magrittr)
 library(parallel)
 #
 # 1.	For each negative probe, calculate total count for all cells. Calculate the median mu.
 neg0 <- test_data$neg0
 mu <- median(rowSums(neg0))

 # 2.	For each positive probe, calculate total count for all cells.
 raw0 <- test_data$raw0
 pos_count <- rowSums(raw0)
 
 # 3.	Select positive probes with total count less than mu, call them low positive probes
 indx_low_pos <- which(pos_count < mu)
 mean(pos_count < mu)
 
 # 4.	Combine low positive probes and negative probes, fit the Poisson Background model, implement the common diagnostics procedure
 negmod <- fitPoisBG(rbind(neg0, raw0[indx_low_pos, ]), size_scale = "sum") # use sum in SMI data. use first will lead to distortion in data
 negdiag2 <- diagPoisBG(negmod, generate_ppplot = FALSE)
 
 
 # 5. perform score tests
 negmod2 <- fitPoisBG(neg0, size_scale = "sum") # use sum in SMI data. use first will lead to distortion in data
 negmod2$sizefact <- negmod$sizefact
 
 sc <- GeoDiff::BGScoreTest(raw0, negmod2, adj = 1, removeoutlier = FALSE, useprior = TRUE)
 # and maybe try different combinations of removeoutlier and useprior
 
 features_high <- ((sc$scores > quantile(sc$scores, probs = 0.4)) & (sc$scores < quantile(sc$scores, probs = 0.95))) %>%
  which() %>%
  names()
 
 # calculate the sizefact
 # estimate a_j from the poisson threshold model
 sizefact0 <- negmod$sizefact
 gamma0 <- mean(negmod2$featfact)
 gamma_features <- rowSums(raw0[features_high,]) - gamma0
 sizefact <- (colSums(raw0[features_high,])-length(features_high)*gamma0*sizefact0)/sum(gamma_features)
 
 # confirm the sum is 1
 sum(sizefact0)
 sum(sizefact)
 
 # percentage of negative
 sum(sizefact<0)/length(sizefact)
 
 # replace negative by 0
 sizefact[sizefact<0] <- 0
 
 # rescale
 gamma_features <- gamma_features*sum(sizefact)
 sizefact <- sizefact/sum(sizefact)
 annot <- test_data$annot
 annot <- as.data.frame(annot)
 rownames(annot) <- colnames(raw0)
 
 annot$fov|>table()
 
 high_ROIs <- names(which(sizefact[annot$fov%in%c(1:2)]>0))
 
 features_all <- rownames(raw0)
 
 NBthDEmod2 <- fitNBthDE(form = ~factor(fov),
                         annot=annot[high_ROIs, ],
                         object=raw0[features_all,high_ROIs],
                         probenum = rep(1, length(features_all)),
                         features_high = features_high,
                         features_all = features_all,
                         sizefact_start=sizefact[high_ROIs],
                         sizefact_BG=sizefact0[high_ROIs],
                         threshold_mean = gamma0,
                         preci2=10000,
                         prior_type="contrast",
                         covrob=FALSE,
                         preci1con=1/25,
                         sizefactrec=FALSE,
                         sizescalebythreshold=TRUE,
                         run_parallel = TRUE,
                         n_parallel = (parallel::detectCores()))
 expect_error(
  NBthDEmod2 <- fitNBthDE(form = ~factor(fov),
                          annot=annot[high_ROIs, ],
                          object=raw0[features_all,high_ROIs],
                          probenum = rep(1, length(features_all)),
                          features_high = features_high,
                          features_all = features_all,
                          sizefact_start=sizefact[high_ROIs],
                          sizefact_BG=sizefact0[high_ROIs],
                          threshold_mean = gamma0,
                          preci2=10000,
                          prior_type="contrast",
                          covrob=FALSE,
                          preci1con=1/25,
                          sizefactrec=FALSE,
                          sizescalebythreshold=TRUE,
                          run_parallel = TRUE,
                          n_parallel = (parallel::detectCores()+1)))
 
})
