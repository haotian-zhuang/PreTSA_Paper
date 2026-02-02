suppressPackageStartupMessages(library(PreTSA))
suppressPackageStartupMessages(library(Matrix))

f <- "sc_pbmc3k"

expr <- readRDS(paste0(f, "/data/expr.rds"))
pseudotime <- readRDS(paste0(f, "/data/pseudotime.rds"))
expr <- expr[, names(pseudotime)]
pseudotime_permute <- readRDS(paste0(f, "/data/pseudotime_permute.rds"))

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res <- temporalTest(expr = expr, pseudotime = pseudotime, pseudotime_permute = pseudotime_permute, knot = 0)
  res[, "gene"] <- rownames(res)
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res, file = paste0(f, "/test/pretsa/res.rds"))
saveRDS(time, file = paste0(f, "/test/pretsa/time.rds"))
saveRDS(memory, file = paste0(f, "/test/pretsa/memory.rds"))
