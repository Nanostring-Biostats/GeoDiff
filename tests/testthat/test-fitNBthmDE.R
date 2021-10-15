### There are two main functions being testing here
### fitNBthmDE with random intercept effect and random slope effect

test_that("fitNBthmDE produces desired results, CTA", {


  #### Specs for fitNBthmDE:
  #   1 The function outputs para.
  #     This matrix has features_all in the columns
  #     and parameters(regression coefficients, threshold, r) in the rows.
  #     Both threshold and r are positive.


  library(Biobase)
  library(dplyr)
  # Preamble/load example data
  data("demoData")
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  # Change to the temporary directory (will set back to getwd() once block is executed.)
  withr::local_dir(tmp_dir)
  expect_true(inherits(demoData, "NanoStringGeoMxSet"))


  # susbet samples
  demoData <- demoData[, c(1:5, 33:37)]
  demoData <- fitPoisBG(demoData, size_scale = "sum")
  demoData <- aggreprobe(demoData, use = "cor")
  demoData <- BGScoreTest(demoData)
  demoData$slidename <- substr(demoData[["slide name"]], 12, 17)
  thmean <- 1 * mean(fData(demoData)$featfact, na.rm = TRUE)
  # Negative and Non-Negative facets:
  demo_pos <- demoData[which(!fData(demoData)$CodeClass == "Negative"), ]
  demo_neg <- demoData[which(fData(demoData)$CodeClass == "Negative"), ]
  # feature factors per groupvar (i.e., slide name)
  sc1_scores <- fData(demo_pos)[, "scores"]
  names(sc1_scores) <- fData(demo_pos)[, "TargetName"]
  # scaling factors
  features_high <- ((sc1_scores > quantile(sc1_scores, probs = 0.4)) &
                      (sc1_scores < quantile(sc1_scores, probs = 0.95))) |>
    which() |>
    names()
  set.seed(123)
  # Negative Binomial threshold model
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
  # fixed effect
  pData(demoData)$group <- c(rep(1, 5), rep(2, 5))
  # run fixed effect model
  NBthDEmod2 <- fitNBthDE(form = ~group,
                          split = FALSE,
                          object = demoData,
                          ROIs_high = ROIs_high,
                          features_high = features_high,
                          features_all = features_all,
                          sizefact_start = demoData[, ROIs_high][['sizefact_fitNBth']],
                          sizefact_BG = demoData[, ROIs_high][['sizefact']],
                          threshold_mean = notes(demoData)[["threshold"]],
                          preci2=10000,
                          prior_type="contrast",
                          covrob=FALSE,
                          preci1con=1/25,
                          sizescalebythreshold=TRUE)


  ### Case 1: run random intercept model
  NBthmDEmod1 <- fitNBthmDE(form = ~ group + (1 | `slide name`),
                            split = FALSE,
                            object = demoData,
                            ROIs_high = ROIs_high,
                            features_all = features_all[1:5],
                            sizefact = demoData[, ROIs_high][["sizefact_fitNBth"]],
                            sizefact_BG = demoData[, ROIs_high][["sizefact"]],
                            preci1=NBthDEmod2$preci1,
                            threshold_mean = thmean,
                            preci2=10000,
                            sizescale = TRUE,
                            controlRandom=list(nu=12, nmh_e=400, thin_e=60))

  # 1: The function outputs para...
  expect_true("para" %in% names(NBthmDEmod1))

  # This matrix has features_all in the columns
  expect_true(all(features_all[1:5] == colnames(NBthmDEmod1$para)))

  # and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "group", "r", "threshold") %in% rownames(NBthmDEmod1$para)))

  # Both threshold and r are positive.
  expect_true(all(NBthmDEmod1$para["threshold",] > 0))
  expect_true(all(NBthmDEmod1$para["r",] > 0))



  ### Case 2: run random slope model
  NBthmDEmod1slope <- fitNBthmDE(form = ~ group + (1 + group | `slide name`),
                            split = FALSE,
                            object = demoData,
                            ROIs_high = ROIs_high,
                            features_all = features_all[1:5],
                            sizefact = demoData[, ROIs_high][["sizefact_fitNBth"]],
                            sizefact_BG = demoData[, ROIs_high][["sizefact"]],
                            preci1=NBthDEmod2$preci1,
                            threshold_mean = thmean,
                            preci2=10000,
                            sizescale = TRUE,
                            controlRandom=list(nu=12, nmh_e=400, thin_e=60))

  # 1: The function outputs para...
  expect_true("para" %in% names(NBthmDEmod1slope))

  # This matrix has features_all in the columns
  expect_true(all(features_all[1:5] == colnames(NBthmDEmod1slope$para)))

  # and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "group", "r", "threshold") %in% rownames(NBthmDEmod1slope$para)))

  # Both threshold and r are positive.
  expect_true(all(NBthmDEmod1slope$para["threshold",] > 0))
  expect_true(all(NBthmDEmod1slope$para["r",] > 0))
})





test_that("fitNBthmDE produces desired results, WTA", {


  ### Same overall workflow as above but with the WTA kidney dataset

  library(dplyr)
  ### Initializing WTA objects before running tests
  # Create temporary directory that will get destroyed after this block is executed.
  tmp_dir <- withr::local_tempdir(pattern = "tmp_dir")
  withr::local_dir(tmp_dir)
  # load kidney data
  data("kidney")
  # subset samples
  kidney <- kidney[, c(1:5, 22:27, 46:50)]
  # alter featureNames
  stopifnot(identical(rownames(fData(kidney)), featureNames(kidney)))
  stopifnot(identical(rownames(fData(kidney)), rownames(exprs(kidney))))
  rownames(fData(kidney))[which(!fData(kidney)$Negative)] <- fData(kidney)[which(!fData(kidney)$Negative), "TargetName"]
  featureNames(kidney) <- rownames(fData(kidney))
  rownames(exprs(kidney)) <- rownames(fData(kidney))

  kidney <- fitPoisBG(kidney, size_scale = "sum")
  kidney <- fitPoisBG(kidney, groupvar = "slide name", size_scale = "sum")
  all0probeidx <- which(rowSums(exprs(kidney))==0)
  kidney <- kidney[-all0probeidx, ]
  kidney <- aggreprobe(kidney, use = "cor")

  # Negative Binomial threshold model
  set.seed(123)
  kidney <- fitNBth(kidney,
                    iterations = 5,
                    start_para = c(200, 1),
                    lower_sizefact = 0,
                    lower_threshold = 100, tol = 1e-8
  )
  # scaling factors
  posdat <- kidney[-which(fData(kidney)$CodeClass == "Negative"), ]
  posdat <- exprs(posdat)
  gene_sum <- rowSums(posdat)
  kidney_neg <- kidney[which(fData(kidney)$CodeClass == "Negative"), ]
  featfact <- fData(kidney_neg)[, "featfact"]
  thmean <- 1*mean(featfact)
  features_high <- ((gene_sum>quantile(gene_sum, probs = 0.5)) & (gene_sum<quantile(gene_sum, probs = 0.95))) |> which() |> names()
  set.seed(123)
  features_high <- sample(features_high, 1500)
  ROIs_high <- sampleNames(kidney)[which((quantile(fData(kidney)[["para"]][, 1],
                                                   probs = 0.90, na.rm = TRUE) - notes(kidney)[["threshold"]])*kidney$sizefact_fitNBth>2)]
  features_all <- rownames(posdat)
  # fixed effect model
  NBthDEmod2 <- fitNBthDE(form = ~region,
                          split = FALSE,
                          object = kidney,
                          ROIs_high = ROIs_high,
                          features_high = features_high,
                          features_all = features_high,
                          sizefact_start = kidney[, ROIs_high][['sizefact_fitNBth']],
                          sizefact_BG = kidney[, ROIs_high][['sizefact']],
                          threshold_mean = notes(kidney)[["threshold"]],
                          preci2=10000,
                          prior_type="contrast",
                          covrob=FALSE,
                          preci1con=1/25,
                          sizescalebythreshold=TRUE)


  ### Case 1: run examplar function from vignette and check each specification
  # random intercept model
  NBthmDEmod2 <- fitNBthmDE(object = kidney,
                            form = ~ region+(1|`slide name`),
                            ROIs_high = ROIs_high,
                            split = FALSE,
                            features_all = features_high[1:5],
                            sizefact = kidney[, ROIs_high][['sizefact_fitNBth']],
                            sizefact_BG=kidney[, ROIs_high][['sizefact']],
                            preci1=NBthDEmod2$preci1,
                            threshold_mean = thmean,
                            preci2=10000,
                            sizescale = TRUE,
                            controlRandom=list(nu=12, nmh_e=400, thin_e=60))

  # 1: The function outputs para...
  expect_true("para" %in% names(NBthmDEmod2))

  # This matrix has features_all in the columns
  expect_true(all(features_high[1:5] == colnames(NBthmDEmod2$para)))

  # and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "regiontubule", "r", "threshold") %in% rownames(NBthmDEmod2$para)))

  # Both threshold and r are positive.
  expect_true(all(NBthmDEmod2$para["threshold",] > 0))
  expect_true(all(NBthmDEmod2$para["r",] > 0))



  ### Case 2: run examplar function from vignette and check each specification
  # random slope model
  NBthmDEmod2slope <- fitNBthmDE(object = kidney,
                                 form = ~ region+(1+region|`slide name`),
                                 split = FALSE,
                                 ROIs_high = ROIs_high,
                                 features_all = features_high[1:5],
                                 sizefact = kidney[, ROIs_high][['sizefact_fitNBth']],
                                 sizefact_BG=kidney[, ROIs_high][['sizefact']],
                                 preci1=NBthDEmod2$preci1,
                                 threshold_mean = thmean,
                                 preci2=10000,
                                 sizescale = TRUE,
                                 controlRandom=list(nu=12, nmh_e=400, thin_e=60))
  # 1: The function outputs para...
  expect_true("para" %in% names(NBthmDEmod2slope))

  # This matrix has features_all in the columns
  expect_true(all(features_high[1:5] == colnames(NBthmDEmod2slope$para)))

  # and parameters(regression coefficients, threshold, r) in the rows.
  expect_true(all(c("(Intercept)", "regiontubule", "r", "threshold") %in% rownames(NBthmDEmod2slope$para)))

  # Both threshold and r are positive.
  expect_true(all(NBthmDEmod2slope$para["threshold",] > 0))
  expect_true(all(NBthmDEmod2slope$para["r",] > 0))
})
