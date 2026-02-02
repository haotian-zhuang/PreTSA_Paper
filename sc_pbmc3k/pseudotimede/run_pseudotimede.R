suppressPackageStartupMessages(library(PseudotimeDE, lib.loc = "/hpc/group/jilab/hz/rlibs"))
suppressPackageStartupMessages(library(Matrix))

f <- "sc_pbmc3k"

expr <- readRDS(paste0(f, "/data/expr.rds"))
pseudotime <- readRDS(paste0(f, "/data/pseudotime.rds"))
expr <- expr[, names(pseudotime)]
pseudotime_permute <- readRDS(paste0(f, "/data/pseudotime_permute.rds"))

ori_tbl <- tibble::tibble(cell = names(pseudotime), pseudotime = unname(pseudotime))
sub_tbl <- lapply(pseudotime_permute, function(i) {
  tibble::tibble(cell = names(i), pseudotime = unname(i))
})

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res <- runPseudotimeDE(gene.vec = rownames(expr), ori.tbl = ori_tbl, sub.tbl = sub_tbl,
                         mat = expr, model = "gaussian", mc.cores = 1)
  res <- as.data.frame(res)
  res[, "gene"] <- rownames(res) <- rownames(expr)
  
  colnames(res)[colnames(res) == "para.pv"] <- "pval"
  colnames(res)[colnames(res) == "emp.pv"] <- "pval.empirical"
  colnames(res)[colnames(res) == "fix.pv"] <- "pval.fixed"
  res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
  res[, "fdr.empirical"] <- stats::p.adjust(res[, "pval.empirical"], method = "fdr")
  res[, "fdr.fixed"] <- stats::p.adjust(res[, "pval.fixed"], method = "fdr")
  res <- res[order(res[, "pval"], res[, "pval.empirical"], res[, "pval.fixed"]), ]
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res, file = paste0(f, "/test/pseudotimede/res.rds"))
saveRDS(time, file = paste0(f, "/test/pseudotimede/time.rds"))
saveRDS(memory, file = paste0(f, "/test/pseudotimede/memory.rds"))
