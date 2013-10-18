# Functions and commands to analyse the Palson data under the hypothesis that gene expression profiles are periodically organised.


# Load the required packages.
library(fMultivar);
library(boot);
library(ggplot2);
library(diptest);
library(gridExtra);


## Generic Functions.

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


## Analysis functions.

neighbour_correlation <- function(expressionMatrix, genePositions, per = FALSE, n = 20, geneSubset = "", cc = c("spearman", "pearson", "kendall")[1], cyclic = TRUE, ...) {
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
  # Get the subset of genes to work with
  if ( nchar(geneSubset) == 0 ) {
    geneList <- rownames(geMat)
  }
  else {
    geneList <- scan(geneSubset, what = "character")
  }
  # Calculating correlations.
  for ( gn in geneList ) {
    geneExPr <- geMat[gn,];
    geneName <- gn;   #TODO change the variable in the rest of the code and test.
    # Checks the existence of the geneName...
    if ( geneName %in% names(geneCoords) ) {
      geneIndex <- which(geneCoords == geneCoords[[geneName]])[1];
    }
    #TODO check this condition! else {stop(sprintf("Can not find gene %s in the positions file", geneName))}
    else { next }
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
  geneCoords <- geneCoords[names(neiCorr)];
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


## Visualisation functions.

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



plot_correlation_densities <- function(posCorr1D, posCorr3D, thres = 0.99, bootData = FALSE, krnl = "epanechnikov", ti = "Position Correlation Densities\n", ...) {
# Plot the density distributions (i.e. histograms of the positional correlations in 1D and 3D.
  # First generate two factors.
  ps1d <- data.frame(posCorr = posCorr1D$corr);
  ps3d <- data.frame(posCorr = posCorr3D$corr);
  ps1d$coords <- "Genomic";
  ps3d$coords <- "SCM";
  # Then combine then in the whole data.
  posCorrs <- rbind(ps1d, ps3d);
  # Do the dip bimodality test and get the p value.
  pv <- dip.test(posCorr3D$corr)$p.value;
  # Plotting part, basic ggplot2
  g <- ggplot(posCorrs, aes(x = posCorr));
  # Densities
  g <- g + geom_density(aes(fill = coords), alpha = 0.75, kernel = krnl, ...);
  # Plot decorations
  g <- g + labs(x = "Position Correlation", y = "Density", title = sprintf("%sDip test p value: %.2e", ti, pv)) + guides(fill = guide_legend(title = "Coord")) + theme(legend.justification = c(1,1), legend.position = c(1,1));
  # Make the title red is pv <= thres
  if ( pv <= (1 - thres) ) {
    g <- g + theme(plot.title = element_text(colour = "red3"));
  }
  # Plot the bootstrap threshold.
  if ( ! identical(bootData, FALSE) ) {
    cutt <- quantile(bootData, probs = thres);
    return(g + geom_vline(xintercept = cutt, col = "red"))
  }
  else {
    return(g)
  }
  # Brilliant ggplot can return a plot object too!!!, if it is not assigned then it will plot it automatically.
}



## Multiscale analysis functions.

## (multiscale analysis entails the analysis of the positional correlation on many different steps that correspond to different sizes of the neighbourhood for correlation (parameter n).

run_bootstrap <- function(exPr, outfileName, nTimes = c(5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100), ...) {
# Run the bootstrap as many times as specified in the "times" list of the number of neighbours.
  bdf <- data.frame();
  for ( i in nTimes ) {
    bs <- bootstrap_data(exPr, i, ...);
    if ( ncol(bdf) != 0 ) {
      bdf[[sprintf("n%s", i)]] <- bs;
    }
    else {
      bdf <- data.frame(bs);
      colnames(bdf) <- sprintf("n%s", i);
    }
  }
  zz <- file(outfileName, "w");
  write.table(bdf, zz);
  close(zz);
}


multi_pos_correlation <- function(exPr, genePos, per, outfileName,  nTimes = c(5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100), ...) {
# Run the position correlation calculations over many different gene neigoubourhoods (as specified by the "nTimes" param.) for a given period.
  # Initialise the data frames.
  posCorDF1D <- data.frame();
  posCorDF3D <- data.frame();

  # Set up a nested function (that does most of the job actually and avoids code duplication).
  posCorr2DF <- function(exPr, genePos, per, n, posCordf, ...) {
  # Internal function to put the positional correlations in a data frame.
    posCor <- neighbour_correlation(exPr, genePos, per, n, ...);
    if ( ncol(posCordf) != 0 ) {
      posCordf[[sprintf("n%s", i)]] <- posCor$corr;
    }
    else {
      # The first column of the data frame will be the gene positions.
      posCordf <- data.frame(posCor$pos);
      colnames(posCordf) <- "pos";
      # Then start populating with the positional correlation score.
      posCordf[[sprintf("n%s", i)]] <- posCor$corr;
    }
    return(posCordf)
  }

  # Loop through nsi
  for ( i in nTimes ) {
    posCorDF1D <- posCorr2DF(exPr, genePos, per = FALSE, n = i, posCordf = posCorDF1D, ...);
    posCorDF3D <- posCorr2DF(exPr, genePos, per = per, n = i, posCordf = posCorDF3D, ...);
  }
  zz1 <- file(sprintf("%s_p%s_1D.out", outfileName, per), "w");
  zz2 <- file(sprintf("%s_p%s_3D.out", outfileName, per), "w");
  write.table(posCorDF1D, zz1);
  write.table(posCorDF3D, zz2);
  close(zz1);
  close(zz2);
}


plot_mulitscale_experiment <- function(posCorr1D_File, posCorr3D_File, bootstrap_File, per, thres = 0.99, ...) {
# Function to create plots out of the results of a multiscale experiment.

  # Nested function to create layouts
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y);

  # Some tricks in the read to make reading faster.
  posCorr1D <- read.table(posCorr1D_File, header = TRUE, sep = " ", comment.char = "", stringsAsFactors = FALSE);
  posCorr3D <- read.table(posCorr3D_File, header = TRUE, sep = " ", comment.char = "", stringsAsFactors = FALSE);
  bootData <- read.table(bootstrap_File, header = TRUE, sep = " ", comment.char = "", stringsAsFactors = FALSE, nrows = 297101, colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"));
  # Reads REALLY fast, however it might be subject to changes, so the old way is kept above.
#  bootData <- fread(bootstrap_File, header = TRUE, autostart = 1)
  # Calculate the best gridding.
  nn <- length(colnames(bootData));
  nc <- ceiling(sqrt(nn));
  if ( nn - floor(sqrt(nn))*ceiling(sqrt(nn)) <= 0 ) {
    nr <- floor(sqrt(nn));
  }
  else {
    nr <- ceiling(sqrt(nn));
  }
  # Set the grid.
  pushViewport(viewport(layout = grid.layout(nr, nc)))
  # bootData is the only data frame where the columns contain ONLY the different neighbouring sweeps. The other two contain a column for the gene positions as well.
  for ( i in 1:length(colnames(bootData)) ) {
    pc1d <- list();
    pc3d <- list();
    n <- colnames(bootData)[i];
    bd <- bootData[[n]];
    pc1d[["corr"]] <- posCorr1D[[n]];
    pc3d[["corr"]] <- posCorr3D[[n]];
    pp <- plot_correlation_densities(pc1d, pc3d, thres, bd, ti = sprintf("Neighbours %s, Period %s\n", n, per), ...);
    # Derive the plot coordinates.
    if ( i%%nc != 0 ) {
      k <- i%/%nc + 1; # The quotient of the division with ncols.
      l <- i%%nc; # The remainder of the division by the ncols.
    }
    else {
      k <- i%/%nc;
      l <- nc;
    }
    print(pp, vp = vplayout(k, l));
  }
}

