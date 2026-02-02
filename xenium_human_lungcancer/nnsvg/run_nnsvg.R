suppressPackageStartupMessages(library(nnSVG))
suppressPackageStartupMessages(library(Matrix))

f <- "xenium_human_lungcancer"

counts <- readRDS(paste0(f, "/data/counts.rds"))
coord <- readRDS(paste0(f, "/data/coord.rds"))
counts <- counts[, rownames(coord)]
expr <- log2(t(t(counts) / colSums(counts) * 1e4) + 1)

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res <- nnSVG(input = expr, spatial_coords = coord, order = "Sum_coords", n_threads = 1)
  res <- as.data.frame(res)
  res <- res[order(res[, "rank"]), ]
  res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
  res[, "gene"] <- rownames(res)
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res, file = paste0(f, "/test/nnsvg/res.rds"))
saveRDS(time, file = paste0(f, "/test/nnsvg/time.rds"))
saveRDS(memory, file = paste0(f, "/test/nnsvg/memory.rds"))
