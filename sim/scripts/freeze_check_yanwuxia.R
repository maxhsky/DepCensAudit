# Precise check of E2 thresholds and the contradiction zone (freeze review)
cat("--- Exact one-sided binomial test H0 p=0.10, alpha=0.05 ---\n")
for (nn in c(2000, 500)) {
  x <- 0:nn
  xmax <- max(x[pbinom(x, nn, 0.10) <= 0.05])
  cat(sprintf("n=%d : exact-test critical x = %d (rate %.4f)\n", nn, xmax, xmax / nn))
}
cat("--- Power at SAP's stated thresholds ---\n")
cat(sprintf("n=2000, x<=173 (0.0865): power at p=0.08 = %.3f\n",
            pbinom(173, 2000, 0.08)))
cat(sprintf("n=2000, x<=177 (exact crit): power at p=0.08 = %.3f\n",
            pbinom(177, 2000, 0.08)))
cat(sprintf("n=500, x<=36 (0.072): power at p=0.06 = %.3f\n",
            pbinom(36, 500, 0.06)))
cat(sprintf("n=500, x<=38 (exact crit): power at p=0.06 = %.3f\n",
            pbinom(38, 500, 0.06)))
cat("--- Contradiction zone: rates passing exact test but failing 0.0866 ---\n")
for (xx in 174:177) {
  cat(sprintf("x=%d (rate %.4f): exact p = %.4f %s 0.05 ; SAP threshold says %s\n",
              xx, xx / 2000, pbinom(xx, 2000, 0.10),
              ifelse(pbinom(xx, 2000, 0.10) <= 0.05, "<=", ">"),
              ifelse(xx / 2000 <= 0.0866, "PASS", "FAIL")))
}
