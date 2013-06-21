# R code to analyse results from the 

periodogram <- function(periodData, pvCuttoff = 0.05, weight = FALSE, tit = "Periods of ...", ...) {
# Plot a human readable period representation of the global spectrum output. (barplot)
  # Sort data according to period.
  d <- periodData[order(periodData[["p"]]),];
  # Get the significant periods.
  if ( weight ) {
    dSelect <- subset(d, d$pv_weighted <= pvCuttoff);
    pv <- dSelect$pv_weighted;
  }
  else {
    dSelect <- subset(d, d$pv <= pvCuttoff);
    pv <- dSelect$pv;
  }
  score <- dSelect$score
  # Calculate a pvalue gradient.
  a <- (score -  min(score)) / (max(score) - min(score));
  # Get the period names.
  nm <- dSelect$p;
  # make the plot.
  barplot(-log10(pv), log = "y", names = nm, col = topo.colors(1, alpha = a), ylim = range(-log10(pv)), xlab = "Period", ylab = "-log(p-value)", main = tit, ...);
  legend("topleft", legend = round(range(score), 2), col = topo.colors(1, alpha = range(a)), pch = 15, pt.cex = 2, horiz = TRUE, title = "Per. Score");
}
