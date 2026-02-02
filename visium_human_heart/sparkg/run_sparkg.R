suppressPackageStartupMessages(library(SPARK))
suppressPackageStartupMessages(library(Matrix))

f <- "visium_human_heart"

counts <- readRDS(paste0(f, "/data/counts.rds"))
coord <- readRDS(paste0(f, "/data/coord.rds"))
counts <- counts[, rownames(coord)]

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res <- sparkg_sw(counts = counts, location = coord)
  res <- res[["res_mtest"]]
  res <- res[order(res[, "combinedPval"]), ]
  res[, "fdr"] <- stats::p.adjust(res[, "combinedPval"], method = "fdr")
  res[, "gene"] <- rownames(res)
  colnames(res)[colnames(res) == "combinedPval"] <- "pval"
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res, file = paste0(f, "/test/sparkg/res.rds"))
saveRDS(time, file = paste0(f, "/test/sparkg/time.rds"))
saveRDS(memory, file = paste0(f, "/test/sparkg/memory.rds"))
