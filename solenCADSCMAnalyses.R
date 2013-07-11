#!/usr/bin/Rscript --vanilla 
# An R function for visualising results from the global spectrum program of the SolenCAD suite.
# Costas Bouyioukos (cbouyio@gmail.com), iSSB, 2013

periodogram <- function(solenCADoutfile, pvThres = 0.05, gTitle = "", pvWeight = FALSE, ...) {
# Plot a human readable periods graph (a barplot or periodogram) from the output file of the global spectrum analysis program.
  # Read data
  perData <- read.table(solenCADoutfile, header = TRUE);
  # Sort data according to period.
  perD <- perData[order(perData[["P"]]),];
  # Get periods that pass the specified p value threshold.
  if ( pvWeight && "pv_weighted" %in% colnames(perD) ) {
    perSelect <- subset(perD, perD$pv_weighted <= pvThres);
    pv <- perSelect$pv_weighted;
    w <- " weighted";
  }
  else if ( pvWeight && !("pv_weighted" %in% colnames(perD)) ) {
    stop("Input Error: Data does not have a weighted p value.");
  }
  else {
    perSelect <- subset(perD, perD$pv <= pvThres);
    pv <- perSelect$pv;
    w <- "";
  }
  # Obtain the score.
  score <- perSelect$score;
  # Calculate a score gradient (to specify the color grading).
  a <- (score -  min(score)) / (max(score) - min(score));
  # Get the period names.
  nm <- perSelect[["P"]];
  # Put together the title.
  gTitle <- sprintf("%s (pv_thres %s%s)", gTitle, pvThres, w);
  # Make the plot.
  barplot(-log10(pv), log = "y", names = nm, col = topo.colors(1, alpha = a), ylim = range(-log10(pv)), xlab = "Period", ylab = "-log(p-value)", main = gTitle, ...);
  # and the legend.
  legend("topleft", legend = round(range(score), 2), col = topo.colors(1, alpha = range(a)), pch = 15, pt.cex = 2, horiz = TRUE, title = "Per. Score");
}

# Set the defaults.
pvThres = 0.05;
gTitle = "";
pvWeight = FALSE;

# When the file is run as a batch.
cat("Global spectrum output file name: ");
outfile <- readLines(con="stdin", 1);

if ( outfile == "" ) {
  cat("Give a global spectrum output file name (no defaults): ");
  outfile <- readLines(con="stdin", 1);
}

cat("Give a p value threshold (default = 0.05): ");
pv <- readLines(con="stdin", 1);
cat("Set an informative title (default = ''): ");
tl <- readLines(con="stdin", 1);

if ( pv != "" ) { pvThres = as.numeric(pv) }
if ( tl != "" ) { gTitle = tl }

quartz();
periodogram(outfile, pvThres, gTitle, pvWeight);
locator();

