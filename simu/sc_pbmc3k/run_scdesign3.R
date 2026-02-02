suppressPackageStartupMessages(library(scDesign3))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

f <- "sc_pbmc3k"
signal <- "v1"

counts <- readRDS(paste0(f, "/data/counts.rds"))
pseudotime <- readRDS(paste0(f, "/data/pseudotime.rds"))

counts <- counts[rowSums(counts > 0) >= 20, ]
set.seed(2025)
counts <- counts[sample(rownames(counts), 2000), names(pseudotime)]

example_sce <- SingleCellExperiment(list(counts = counts),
                                    colData = data.frame(row.names = names(pseudotime),
                                                         pseudotime = pseudotime))

set.seed(1)
example_data <- construct_data(
  sce = example_sce,
  assay_use = "counts",
  celltype = NULL,
  pseudotime = "pseudotime",
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "ind"
)

example_marginal <- fit_marginal(
  data = example_data,
  predictor = "gene",
  mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
  sigma_formula = "s(pseudotime, k = 5, bs = 'cr')",
  family_use = "nb",
  n_cores = 20,
  usebam = FALSE
)

set.seed(1)
example_copula <- fit_copula(
  sce = example_sce,
  assay_use = "counts",
  marginal_list = example_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 20,
  input_data = example_data$dat
)

example_para <- extract_para(
  sce = example_sce,
  marginal_list = example_marginal,
  n_cores = 20,
  family_use = "nb",
  new_covariate = example_data$newCovariate,
  data = example_data$dat
)

diff <- apply(example_para$mean_mat, 2, function(x){max(x)-min(x)})
diff_ordered <- order(diff, decreasing = TRUE)
num_de <- round(length(diff_ordered)/10)
ordered <- diff[diff_ordered]
if (signal == "v1") {
  de_idx <- names(ordered)[1:num_de] # top 10%
  non_de_idx <- names(ordered)[-(1:num_de)]
} else if (signal == "v2") {
  de_idx <- names(ordered)[(1:num_de) + 2*num_de] # top 20% - 30%
  non_de_idx <- names(ordered)[-((1:num_de) + 2*num_de)]
} else if (signal == "v3") {
  de_idx <- names(ordered)[(1:num_de) + 4*num_de] # top 40% - 50%
  non_de_idx <- names(ordered)[-((1:num_de) + 4*num_de)]
} else if (signal == "v4") {
  de_idx <- names(ordered)[(1:num_de) + 6*num_de] # top 60% - 70%
  non_de_idx <- names(ordered)[-((1:num_de) + 6*num_de)]
} else if (signal == "v5") {
  de_idx <- names(ordered)[(1:num_de) + 8*num_de] # top 80% - 90%
  non_de_idx <- names(ordered)[-((1:num_de) + 8*num_de)]
}
###
de_label <- c(rep("tvg", length(de_idx)), rep("nontvg", length(non_de_idx)))
names(de_label) <- c(de_idx, non_de_idx)
saveRDS(de_label, file = paste0("simu/", f, "/simdata/genelist.rds"))
###
non_de_mat <- apply(example_para$mean_mat[,non_de_idx], 2, function(x){
  avg <- mean(x)
  new_mean <- rep(avg, length(x))
  return(new_mean)
})
example_para$mean_mat[,non_de_idx] <- non_de_mat

set.seed(1)
example_newcount <- simu_new(
  sce = example_sce,
  mean_mat = example_para$mean_mat,
  sigma_mat = example_para$sigma_mat,
  zero_mat = example_para$zero_mat,
  quantile_mat = NULL,
  copula_list = example_copula$copula_list,
  n_cores = 20,
  family_use = "nb",
  input_data = example_data$dat,
  new_covariate = example_data$newCovariate,
  important_feature = rep(TRUE, dim(example_sce)[1]),
  filtered_gene = example_data$filtered_gene
)

simu_sce <- SingleCellExperiment(list(counts = example_newcount), colData = example_data$newCovariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

sim_counts <- counts(simu_sce)
saveRDS(sim_counts, file = paste0("simu/", f, "/simdata/sim_counts.rds"))
