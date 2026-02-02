suppressPackageStartupMessages(library(SPARK))
suppressPackageStartupMessages(library(Matrix))

f <- "visium_human_heart"

counts <- readRDS(paste0(f, "/data/counts.rds"))
coord <- readRDS(paste0(f, "/data/coord.rds"))
counts <- counts[, rownames(coord)]

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  obj <- CreateSPARKObject(counts = counts, location = as.data.frame(coord), percentage = 0, min_total_counts = 0)
  obj@lib_size <- colSums(obj@counts)
  obj <- spark.vc(object = obj, covariates = NULL, lib_size = obj@lib_size, fit.model = "poisson", num_core = 1)
  res <- spark.test(object = obj)
  res <- res@res_mtest
  res <- res[order(res[, "combined_pvalue"]), ]
  res[, "fdr"] <- stats::p.adjust(res[, "combined_pvalue"], method = "fdr")
  res[, "gene"] <- rownames(res)
  colnames(res)[colnames(res) == "combined_pvalue"] <- "pval"
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res, file = paste0(f, "/test/spark/res.rds"))
saveRDS(time, file = paste0(f, "/test/spark/time.rds"))
saveRDS(memory, file = paste0(f, "/test/spark/memory.rds"))
