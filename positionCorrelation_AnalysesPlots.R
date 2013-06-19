# Functions and commands to analyse the Palson data under the hypothesis that gene expression profiles are periodically organised.


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
    stop("For the 3D correlation a period must be specified.")
  }
  modCoords <- pos - pos%/%per * per;
  return(modCoords)
}



neighbour_correlation <- function(expressionMatrix, genePositions, n = 10, cc = c("spearman", "pearson", "kendall")[1], per = FALSE, geneList = FALSE, cyclic = TRUE, ...) {
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
      geneIndex <- which(geneCoords == geneCoords[[geneName]])[1]; # TODO fix this with the use of match.
    }
    else {next}
    # Check the positioning of the gene index and generate the neighbourhood.
    if ( geneIndex + n <= noGenes & geneIndex - n >= 0 ) {
      neighCoords <- c(seq(geneIndex - n, length = n), seq(geneIndex + 1, length = n));
    }
    else if ( geneIndex + n > noGenes & geneIndex - n >= 0 ) {
      neighCoords <- c(seq(geneIndex - n, length = n), seq(geneIndex + 1, length = noGenes - geneIndex), seq(1, length = n - (noGenes - geneIndex)));
    }
    else if ( geneIndex + n <= noGenes & geneIndex - n < 0 ) {
      neighCoords <- c(seq(1, length = geneIndex - 1), seq(noGenes - (n - geneIndex + 1), length = n - geneIndex + 1), seq(geneIndex + 1, length = n));
    }
    avgExprProf <- colMeans(geMat[neighCoords,]);
    corr <- cor(colMeans(geneExPr), avgExprProf, method = cc, ...);
    neiCorr <- append(neiCorr, corr);
    geneNames <- append(geneNames, geneName);
  }
  names(neiCorr) <- geneNames;
  neiCorr <- neiCorr[names(geneCoords)];
  # Whather a specific list of genes is asked to be ploted.
  if ( geneList != FALSE ) {
    geneL <- read.table(geneList)[[1]];
    geneNamesL <- intersect(geneNames, geneL);
    geneCoords <- geneCoords[geneNamesL];
    neiCorr <- neiCorr[names(geneCoords)];
  }
  return(list(corr = neiCorr, pos = geneCoords))
}



bootstrap_data <- function(expressionMatrix, bootTimes = 100, n = 10, pv = 0.02, ...) {
# Do bootstraping of the data to obtain a significant spearman correlation threshold.
  noGenes <- dim(expressionMatrix)[1];
  corTld  <- NULL;
  for ( i in 1:bootTimes ) {
    randomOrder <- sample(1:noGenes, noGenes, replace = FALSE);
    newData <- expressionMatrix[randomOrder, ];
    corTot <- cor(t(newData), method = "spearman");
    corMatrix <- matrix(0, ncol = 2*n + 1, nrow = noGenes - 2*n);
    for( i in 1:(noGenes - 2*n) ) {
      corMatrix[i, ] <- corTot[i + n, i:(i + 2*n)];
    }
    corTld <- cbind(corTld, rowSums(corMatrix[, -(n + 1)]) - 1);
  }
  threshold <- quantile(as.vector(corTld), probs = pv);
  return(threshold)
}

