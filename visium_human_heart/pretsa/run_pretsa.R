suppressPackageStartupMessages(library(PreTSA))
suppressPackageStartupMessages(library(Matrix))

f <- "visium_human_heart"

expr <- readRDS(paste0(f, "/data/expr.rds"))
coord <- readRDS(paste0(f, "/data/coord.rds"))
expr <- expr[, rownames(coord)]

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res <- spatialTest(expr = expr, coord = coord, knot = 0)
  res[, "gene"] <- rownames(res)
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res, file = paste0(f, "/test/pretsa/res.rds"))
saveRDS(time, file = paste0(f, "/test/pretsa/time.rds"))
saveRDS(memory, file = paste0(f, "/test/pretsa/memory.rds"))
