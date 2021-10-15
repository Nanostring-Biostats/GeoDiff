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

