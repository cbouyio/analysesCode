# Functions and commands to analyse the Palson data under the hypothesis that gene expression profiles are periodically organised.


# Load the required packages.
library(fMultivar);
library(boot);
library(ggplot2);
library(diptest);

## Functions.

find_gene_name <- function(lst, affProb, ...) {
# Return the common gene name of an Affymetrix probe name.
  nm <- vector();
  for ( i in 1:length(lst) ) {
    if ( ! is.na(lst[[i]][1]) ) {
      if ( lst[[i]][1] == affProb ) {
        nm <- append(nm, names(lst)[[i]]);
      }
    }
  }
  return(nm)
}



modulo_coordinates <- function(pos, per, ...) {
# Return the modulo coordinates of a list of coordinates.
  if ( per == FALSE ) {
    stop("For the solenoid coordinates a period must be specified.")
  }
  modCoords <- pos - pos%/%per * per;
  return(modCoords)
}



neighbour_correlation <- function(expressionMatrix, genePositions, per = FALSE, n = 10, cc = c("spearman", "pearson", "kendall")[1], cyclic = TRUE, ...) {
# Calculate the correlation coefficient between the gene expression profiles of a gene and the average of its N neighbours.
  neiCorr <- vector();
  geneNames <- vector();
  # Read files and generate the data structures.
  geMat <- read.table(expressionMatrix, header = TRUE);
  gp <- read.table(genePositions);
  gpp <- gp[[2]];
  names(gpp) <- gp[[1]];
  # Keep only the genes for which we have expression profiles from the gene expression matrix.
  genePos <- sort(gpp[rownames(geMat)]);
  noGenes <- length(genePos);
  # Transform the position to circular coordinates (if is asked).
  if ( per ) {
    geneCoords <- sort(modulo_coordinates(genePos, per));
  }
  else {
    geneCoords <- genePos;
  }
  # Calculating correlations.
  for ( i in 1:nrow(geMat) ) {
    geneExPr <- geMat[i,];
    geneName <- rownames(geneExPr)[1];
    # Checks the existence of the geneName...
    if ( geneName %in% names(geneCoords) ) {
      geneIndex <- which(geneCoords == geneCoords[[geneName]])[1];
    }
    else {next}
    # Check the positioning of the gene index and generate the neighbourhood.
    if ( geneIndex + n <= noGenes & geneIndex - n > 0 ) {
      neighCoords <- c(seq(geneIndex - n, length = n), seq(geneIndex + 1, length = n));
    }
    else if ( geneIndex + n > noGenes & geneIndex - n > 0 ) {
      neighCoords <- c(seq(geneIndex - n, length = n), seq(geneIndex + 1, length = noGenes - geneIndex), seq(1, length = n - (noGenes - geneIndex)));
    }
    else if ( geneIndex + n < noGenes & geneIndex - n <= 0 ) {
      neighCoords <- c(seq(1, length = geneIndex - 1), seq(noGenes - (n - geneIndex), length = n - geneIndex + 1), seq(geneIndex + 1, length = n));
    }
    # Collect the names of the neighbouring genes.
    neiNames <- names(geneCoords[neighCoords]);
    # Get the gene expression profiles of neighbouring genes.
    neiExprM <- geMat[neiNames,];
    corr <- sum(rapply(as.list(data.frame(t(neiExprM))), cor, x=as.double(geneExPr), method = "spearman"));
    neiCorr <- append(neiCorr, corr);
    geneNames <- append(geneNames, geneName);
  }
  # Assigne the names.
  names(neiCorr) <- geneNames;
  neiCorr <- neiCorr[names(geneCoords)];
  return(list(corr = neiCorr, pos = geneCoords))
}



bootstrap_data <- function(expressionProfiles, n, bootTimes = 100, ...) {
# Do bootstraping of the data to obtain a significant spearman correlation threshold.
  geMat <- read.table(expressionProfiles, header = TRUE);
  noGenes <- nrow(geMat);
  corBoots  <- rep(NA, bootTimes * noGenes);
  k <- 1;
  for ( i in 1:bootTimes ) {
    for ( j in 1:noGenes ) {
      geneExpr <- as.double(geMat[j,]);
      randomSample <- sample(c(1:noGenes)[!c(1:noGenes) == j], n, replace = FALSE);
      rndExprM <- geMat[randomSample,];
      corr <- sum(rapply(as.list(data.frame(t(rndExprM))), cor, x = geneExpr, method = "spearman"))
      corBoots[k] <- corr;
      k <- k + 1;
    }
  }
  return(corBoots)
}



plot_position_correlation <- function(posCorr1D, posCorr3D, btstrp, per, n = "30", pv = 0.001, bins = 35, ...) {
  # Specify the limits of the y axis.
  yl <- range(c(posCorr1D$corr, posCorr3D$corr));
  # Calculate the threshold for the given p-value.
  thres <- quantile(btstrp, probs = 1 - pv);
  # Do the plotting, first the 1D scatterplot.
  par(mfrow = c(2,2), cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8);
  plot(posCorr1D$pos, posCorr1D$corr, xlab = "Genome Coordinate", ylab = "Neighbour Correlation", main = sprintf("Position correlations scatterplot (all genes, n = %s)", n), pch = ".", ylim = yl);
  abline(h = thres);
#  quartz();
  # The 3D scatterplot.
  plot(posCorr3D$pos, posCorr3D$corr, xlab = sprintf("Modulo Coordinate (per = %s)", per), ylab = "Neighbour Correlation", main = sprintf("Position correlations scatterplot (all genes, n = %s)", n), pch = ".", ylim = yl);
  abline(h = thres);
#  quartz();
  # The 1D binned scatterplot.
  plot(squareBinning(posCorr1D$pos, posCorr1D$corr, bin = bins), addRug=FALSE, col = topo.colors(1, alpha = seq(0, 1, length = 24)), xlab = "Genome Coordinate", ylab = "Neighbour Correlation", main = sprintf("Position correlations bin-plot (all genes, n = %s)", n), ylim = yl);
  abline(h = thres);
#  quartz();
  # The 3D binned scatterplot.
  plot(squareBinning(posCorr3D$pos, posCorr3D$corr, bin = bins), addRug=FALSE, col = topo.colors(1, alpha = seq(0, 1, length = 24)), xlab = sprintf("Modulo Coordinate (per = %s)", per), ylab = "Neighbour Correlation", main = sprintf("Position correlation bin-plot (all genes, n = %s)", n), ylim = yl);
  abline(h = thres);
}



plot_correlation_densities <- function(posCorr1D, posCorr3D, thres = 0.95, ...) {
# Plot the density distributions (i.e. histograms of the positional correlations in 1D and 3D.
  # First generate two factors.
  ps1d <- data.frame(posCorr = posCorr1D$corr);
  ps3d <- data.frame(posCorr = posCorr3D$corr);
  ps1d$coords <- "1D";
  ps3d$coords <- "3D";
  # Then combine then in the whole data.
  posCorrs <- rbind(ps1d, ps3d);
  # Do the plot.
  pv <- dip.test(posCorr3D$corr)$p.value;
  geom_density(kernel = "rectangular");
  qplot(posCorr, data = posCorrs, geom = "density", fill = coords, alpha = I(0.75), xlab = "Position Correlation", ylab = "Density", main = sprintf("Position Correlation Distributions\nDip test p value: %.2e", pv));
#  g <- ggplot(posCorr, data = posCorrs, aes(x = posCorr), fill = coords, alpha = I(0.75), xlab = "Position Correlation", ylab = "Density", main = sprintf("Position Correlation Distributions\nDip test p value: %.2e", pv));
#  g + geom_density(kernel = "rectangular");
## Attempts to plot with a rectangular kernel.
}

