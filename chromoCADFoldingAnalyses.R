# 3D plotting and statistics for the folding program output.

library("rgl");
library("stringr");
library("fpc");
library("plyr");
library("MASS");
library("RColorBrewer");
library("ggplot2");
library("gridExtra");
library("lattice");
library("scatterplot3d");



###############################################################################
## Low level functions - Access, store and basic calculations on folding data.
###############################################################################

# Read data from singular files as well as collections (directories etc.).

getCoordinatesData <- function(fileName, point = c("middle", "end")[1], type = c("circular", "linear")[1]) {
#TODO allow linear chromosomes.
# Return data as a frame with the 3D coordinates of the end-points of each cylinder.
# The cylinder number becomes the row name in the DF, and will be used to identify the point.
  d <- read.table(fileName, header = TRUE, comment.char = "&");
  rn <- as.character(d$X.mono);
  # Correction for circular molecules for the last cylinder which has the same coordinates with the first.
  if (rn[1] == rn[length(rn)]) {
    rn[length(rn)] <- paste(rn[length(rn)], "-c", sep = "");
  }
  if ( point == "middle" ) {
    coordsDF <- data.frame(d[2:4], row.names = rn);
  }
  else if ( point == "end" ) {
    coordsDF <- data.frame(d[5:7], row.names = rn);
  }
  return(coordsDF)
}


getMultipleData <- function(dirList, funct, type) {
# Return a list of calculated data from the OUT (the output simulation) directories.
# funct : Specifies the function according to which the measure will be calculated (e.g. "girth")
# type  : Specifies the data upon which the measure will be calculated.
  dfm <- numeric();
  calculateData <- NULL;
  # Get the right function
  if (funct == "pmd") {
    calculateData <- match.fun(meanPairwiseDist);
  }
  else if (funct == "girth") {
    calculateData <- match.fun(gyrationRadius);
  }
  # Get the right data dir
  if (type == "binding") {
    suff <- "Binding_sites_#1_ch_0";
  }
  else if (type == "coordinates") {
    suff <- "Coordonnees_ADN_ch_0";
  }
  # Get the calculated data.
  for (dr in dirList) {
    directory <- sprintf("%s/%s", dr, suff);
    fl <- dir(directory, full.names = TRUE)[1];
    dfm <- append(dfm, calculateData(fl));
  }
  return(dfm)
}


getTracesData <- function(directory, funct) {
# Multivalent function, returns a data frame of the desired collection of simulation files with data calculated by the supported functions of "funct"
# Operates on RES directories.
# Return something equivalent to a time series.
  mcs           <- numeric();
  dfm            <- numeric();
  calculateData <- NULL;
  if (funct == "pmd") {
    calculateData <- match.fun(meanPairwiseDist);
  }
  else if (funct == "girth") {
    calculateData <- match.fun(gyrationRadius);
  }
  else if (funct == "clust") {
    calculateData <- match.fun(clusterDistance);
  }
  listOfFilenames <- dir(sprintf("%s", directory), full.names = TRUE);
  for (f in listOfFilenames) {
    mcs <- append(mcs, as.numeric(str_extract(f, "[0-9]{11}")));
    dfm <- append(dfm, calculateData(f));
  }
  # Return two different data frames depending on the function (actually only the name of the DF component changes).
  if (funct == "pmd") {
    return(list(mcSteps = mcs, pmd = dfm))
  }
  else if (funct == "girth") {
    return(list(mcSteps = mcs, girth = dfm))
  }
  else if (funct == "clust") {
    return(list(mcSteps = mcs, clust = dfm))
  }
}


getMultipleTraceData <- function(directoryList, funct, type, familyIndex = 1) {
# Return a list of traces from multiple data directories, as well as the min and max values of this quantity.
# funct       : Specifies the function according to which the measure will be calculated (e.g. "girth")
# type        : Specifies the data upon which the measure will be calculated.
# familyIndex : The number of different interaction sites families. (i.e. the number of different Binding_sites_#*_ch_0 directories)
  listData <- list();
  maxQuant <- 0;
  minQuant <- Inf;
  if (type == "binding") {
    suff <- sprintf("Binding_sites_#%s_ch_0", str(familyIndex));
  }
  else if (type == "coordinates") {
    suff <- "Coordonnees_ADN_ch_0";
  }
  for (i in 1:length(directoryList)) {
    directory <- sprintf("%s%s", directoryList[i], suff);
    dfm <- getTracesData(directory, funct);
    listData[[i]] <- dfm;
    maxQuant <- max(maxQuant, max(dfm[[funct]]));
    minQuant <- min(minQuant, min(dfm[[funct]]));
  }
  rData <- list(data = listData, min = minQuant, max = maxQuant);
  return(rData)
}


getRatios <- function(dataCollectionPMD, dataCollectionGirth) {
# Calculate a list of ratios between PMDs and Girths
  ratios <- list();
  minRatio <- Inf;
  maxRatio <- 0;
  for (i in 1:length(dataCollectionPMD[["data"]])) {
    ratios[[i]] <- dataCollectionPMD[["data"]][[i]]$pmd / dataCollectionGirth[["data"]][[i]]$girth[-2:0]; # Check the -2 subscripting because the girth calculations include two extraneous files of the 0 and the 1st MC step.
    minRatio <- min(minRatio, min(ratios[[i]]));
    maxRatio <- max(maxRatio, max(ratios[[i]]));
  }
  return(list(ratioList = ratios, min = minRatio, max = maxRatio))
}


# Low level clustering functions.

dbscanClustering <- function(dframe, e = 90, cp = 2, ...) {
# Return named cluster vector showing the cluster assignment of each data point.
# Row names in a properly read dframe should correspond to the number of the cylinder that an interaction site is positioned.
  dbcl <- dbscan(dframe, eps = e, MinPts = cp, ...);
  clusterVector <- dbcl$cluster;
  centres <- double();
  labels <- list();
  for (i in unique(clusterVector)) {
    clusterCentre <- colMeans(dframe[clusterVector == i, ]);
    centres <- rbind(centres, clusterCentre);
    # TODO the labels list....
  }
  names(clusterVector) <- rownames(dframe);
  centresWeighted <- centres %*% abs(princomp(scale(dframe))$loadings[,1]);
  return(list(clusters = clusterVector, centres = centres, weighted.centres = centresWeighted))
}


withinClusterSS <- function(dframe, ... ) {
# Return the intra-cluster sum of square distances out of a data frame of cluster points.
  centroid <- colMeans(dframe);            # The centroid of the cluster
  dtCent <- rbind(dframe, centroid);       # A DF with the points as well as the centroid coordinates
  dists <- as.matrix(dist(dtCent));        # The distance matrix of the above DF
  clustDist <- sum(dists[nrow(dists), ] ** 2); # The mean of the last line of the distance matrix is the mean distance of the cluster points from the centroid
  return(clustDist)
}


clusterDistanceFamilies <- function(dirName, fileIndex = 1) {
# fileIndex : An integer specifying the file that needs to be retrieved from each sorted file list
#             of a directory.
# Return a vector of all the mean euclidean distances of cluster points from the centroid for clusters of the same family.
  clustDistList <- double();
  for (dr in dir(path = dirName, pattern = "Binding_sites#*", full.names = TRUE)) {
    fileName <- dir(dr, full.names = TRUE)[fileIndex];
    dfm <- getCoordinatesData(fileName);
    clustDist <- clusterDistanceDF(dfm);
    clustDistList <- append(clustDistList, clustDist);
    }
  return(clustDistList);
}


clustDistDiff <- function(dirName1, dirName2, fileIndex = 1) {
# Return the distance between two different cluster distance vectors.
# (This measure is a proxy of the degree of "clustering" between two different "clusterings")
  cl1 <- clusterDistance(dirName1, fileIndex);
  cl2 <- clusterDistance(dirName2, fileIndex);
  return(dist(rbind(cl1, cl2)))
}


intraClusterSS <- function(dframe, e = 90, cp = 2, ...) {
# Return a vector (of length equal to the number of clusters) of the intra-cluster sum of squares of each cluster.
  clusterVector <- dbscanClustering(dframe, e, cp, ...)$clusters;
  noClusters <- max(range(clusterVector));
  clusterDists <- vector();
  for (i in seq(noClusters)) {
    dfCluster <- dframe[clusterVector == i, ]; # Get a subset of the data frame with the points corresponding to the cluster i.
    icSS <- withinClusterSS(dfCluster)/nrow(dfCluster);
    clusterDists <- append(clusterDists, icSS);
  }
  names(clusterDists) <- sprintf("Cluster-%i", seq(noClusters));
  return(clusterDists)
}


intraClusterDensity <- function(dframe, e = 90, cp = 2, ...) {
# Return a vector (of length equal to the number of clusters) of the density of each cluster.
  clusterVector <- dbscanClustering(dframe, e, cp, ...)$clusters;
  noClusters <- max(range(clusterVector));
  clusterDens <- vector();
  for (i in seq(noClusters)) {
    dfCluster <- dframe[clusterVector == i, ]; # Get a subset of the data frame with the points corresponding to the cluster i.
    icd <- cloudPointDensity(dfCluster);
    clusterDens <- append(clusterDens, icd);
  }
  names(clusterDens) <- sprintf("Cluster-%i", seq(noClusters));
  return(clusterDens)
}


intraClusterMPD <- function(dframe, e = 90, cp = 2, ...) {
# Return a vector (of length equal to the number of clusters) of the MPD of each cluster.
  clusterVector <- dbscanClustering(dframe, e, cp, ...)$clusters;
  noClusters <- max(range(clusterVector));
  clusterMPD <- vector();
  for (i in seq(noClusters)) {
    dfCluster <- dframe[clusterVector == i, ]; # Get a subset of the data frame with the points corresponding to the cluster i.
    icd <- mean(dist(dfCluster));
    clusterMPD <- append(clusterMPD, icd);
  }
  names(clusterMPD) <- sprintf("Cluster-%i", seq(noClusters));
  return(clusterMPD)
}
# TODO actually the three function above is a code duplication and they should be merged to one generic function with two different function callers.


# General purpose functions.

sphericalVolume <- function(dframe, ...) {
# Return the volume of a sphere with diameter the distance between the two more distant point in the data frame.
  maxDistance <- max(dist(dframe));
  return( (pi*((maxDistance / 1000)^3)) / 6.0)
}


cloudPointDensity <- function(dframe, ...) {
# Return a proxy measure of the density of a cloud of points.
  v <- sphericalVolume(dframe);
  p <- nrow(dframe);
  return(p / v)
}


jaccardIndex <- function(labelsA, labelsB, ...) {
# Return the Jaccard similarity index between two sets.
  u <- union(labelsA, labelsB);
  i <- intersect(labelsA, labelsB);
  return(i/u) # The Jaccard similarity index
}


meanPairwiseDist <- function(fileName) {
# Return the mean pairwise distance between all the points in a single binding file.
  df <- getCoordinatesData(fileName);
  meanDist <- mean(dist(df));
  return(meanDist)
}


gyrationRadius <- function(fileName) {
# Return the radius of gyration (formula taken -but adapted- from wikipedia article on radius of gyration).
  df <- getCoordinatesData(fileName);
  n <- nrow(df);
  sumSq <- sum(dist(df)**2);
  gyration <- sqrt(sumSq/((n**2)));
  return(gyration)
}



###############################################################################
## High level functions - Analyse and plot chromosome folding single file data.
###############################################################################

plotFibre <- function(fileName, lineWidth = 3, subtitle = "Demo", op = FALSE, ...) {
# Baseline 3d plotting function of the fibre.
  if (op == TRUE) {
    open3d(windowRect = c(windowRect = c(0, 0, 800, 800)));
  }
  df <- getCoordinatesData(fileName, point = "end");
  mcStep <- as.numeric(str_extract(fileName, "[0-9]{11}"));
  plot3d(df, type = "l", lwd = lineWidth, xlab = "x or 1", ylab = "y or 2", zlab = "z or 3", box = FALSE, top = TRUE, col = "darkgray", ...);
  title3d(main = sprintf("MC Step: %s", mcStep), sub = subtitle);
}


plotFibreInteractions <- function(directoryName, sizeLine = 2, sizePoint = 5, ...) {
# Plot the fibre as well as the interaction sites.
# Works only for interaction site families on one chromosome!
  tlt <- substr(directoryName, 1, nchar(directoryName) - 1);
  tlt <- tail(strsplit(tl, "_")[[1]], n = 1);
  bindingSitesDirsList <- Sys.glob(sprintf("%s/Binding_sites_#*", directoryName));
  fileNameFiber <- dir(sprintf("%s/Coordonnees_ADN_ch_0/", directoryName), full.names = TRUE);
  plotFibre(fileNameFiber, sizeLine, subtitle = "", main =  sprintf("Fibre-interaction sites plot %s", tlt), ...);
  # Get the number of interacting families to assign different colour.
  cols <- rainbow(length(bindingSitesDirsList));
  for (i in 1:length(bindingSitesDirsList)) {
    fileNameInter <- dir(sprintf("%s/Binding_sites_#%i_ch_0/", directoryName, i), full.names = TRUE);
    dfInter <- getCoordinatesData(fileNameInter);
    points3d(dfInter$rfin.0., dfInter$rfin.1., dfInter$rfin.2., size = sizePoint, col = cols[i], ...);
    texts3d(dfInter$rfin.0., dfInter$rfin.1., dfInter$rfin.2., texts = rownames(dfInter), adj = c(-0.1,-0.1), ...);
  }
}


plotFibreReport <- function(dirName, lineSize = 2, ...) {
# A 2-dimensional 3D plot. (just for the reporting.
  # some code duplication with the previous function.
  tlt <- substr(dirName, 1, nchar(dirName) - 1);
  tlt <- tail(strsplit(tlt, "_")[[1]], n = 1);
  bindingSitesDirsList <- Sys.glob(sprintf("%s/Binding_sites_#*", dirName));
  fileNameFiber <- dir(sprintf("%s/Coordonnees_ADN_ch_0/", dirName), full.names = TRUE);
  dFibre <- getCoordinatesData(fileNameFiber, point = "end");
  sc <- scatterplot3d(dFibre, type = "l", box = FALSE, lwd = lineSize, xlab = "x-axis", ylab = "y-axis", zlab = "z-axis", main = paste("Fibre-interaction sites plot", tlt, sep = " "), lty.grid = "dotted", color = "darkgrey", ...);
  # Different colours for different interacting families.
  cols <- rainbow(length(bindingSitesDirsList));
  for (i in 1:length(bindingSitesDirsList)) {
    fileNameInter <- dir(sprintf("%s/Binding_sites_#%i_ch_0/", dirName, i), full.names = TRUE);
    dInter <- getCoordinatesData(fileNameInter);
    sc$points(dInter, pch = 16, col = cols[i], ...);
    txt2d <- sc$xyz.convert(dInter);
    text(txt2d$x, txt2d$y, labels = rownames(dInter), cex = 0.75, pos = 3);
  }
}


plotFibreBoxplots <- function(outDirListPer, outDirListRnd, funct, type, colours = c("turquoise", "tomato"), ...) {
# Returns two box-plots from two collections of OUT files.
  if (funct == "pmd") {
    label <- "Mean Pairwise Distance"
  }
  else if (funct == "girth") {
    label <- "Gyration Radius"
  }
  else {
    stop("Function statement unsupported, provide one of 'pmd' or 'girth'");
  }
  perData <- getMultipleData(outDirListPer, funct, type);
  rndData <- getMultipleData(outDirListRnd, funct, type);
  boxplot(perData, rndData, varwidth = TRUE, names = c("Periodic", "Random"), ylab = "Nanometres", main = sprintf("%s Boxplots", label), boxwex = 0.6, col = colours, ...);
  legend("bottomright", legend = sprintf("U-test: %.2f\np value: %.4f", wilcox.test(perData, rndData)$statistic, wilcox.test(perData, rndData)$p.value))
}


plotFibreTraces <- function(directoryName, funct, ...) {
# Plots the simulation traces of either the girth or the pmd, depending on the "funct" call.
  if (funct == "pmd")  {
    label <- "Mean Pairwise Distance"
  }
  else if (funct == "girth") {
    label <- "Gyration Radius"
  }
  else {
    stop("Function statement unsupported, provide one of 'pmd' or 'girth'");
  } 
  dfm <- getTracesData(directoryName, funct);
  plot(dfm$mcSteps, dfm[[funct]], type = 'l', xlab = "MC step", ylab = label, ...);
}


# Functions relative to the clustering analysis.

clusterNumber <-function(dirList, ...) {
# Return a vector of the lenght of the dirList with the number of clusters of each experiment.
  nclusters <- vector();
  for ( dr in dirList ) {
    directory <- sprintf("%s/Binding_sites_#1_ch_0", dr);
    fileName <- dir(directory, full.names = TRUE)[1];
    dfm <- getCoordinatesData(fileName);
    nc <- max(dbscanClustering(dfm)$clusters);
    nclusters <- append(nclusters, nc);
  }
  return(nclusters)
}


plotClusterNumbers <- function(dirListPer, dirListRnd, type = c("bar", "box")[1], jitt = FALSE, ...) {
# Plot the barplot of cluster numbers plus error bars.
  ncPer <- clusterNumber(dirListPer);
  ncRnd <- clusterNumber(dirListRnd);
  lenP <- length(ncPer);
  lenR <- length(ncRnd);
  dfm <- data.frame(Arrangement = factor(c(rep("Periodic", lenP), rep("Random", lenR))), NoClusters = c(ncPer, ncRnd));
  uTest <- wilcox.test(ncPer, ncRnd);
  p1 <- ggplot(dfm, aes(fill = Arrangement, factor(NoClusters)));
  p1 <- p1 + geom_bar(position = "dodge") + labs(x = "Number of Clusters", y = "Count") ;
  p2 <- ggplot(dfm, aes(Arrangement, NoClusters));
  p2 <- p2 + geom_boxplot(aes(fill = Arrangement))  + labs(y = "Number of Clusters");
  if ( jitt == TRUE) {
    p2 <- p2 + geom_jitter();
  }
  if ( type == "bar" ) {
    return(p1 + ggtitle(sprintf("U-Test p value: %.3e", uTest$p.value)));
  }
  else if ( type == "box" ) {
    return(p2 + ggtitle(sprintf("U-Test p value: %.3e", uTest$p.value)));
  }
}


getIntraCusterMeasure <- function(dirList, funct, ...) {
# Return a vector of the required intra-cluster measure for many directories.
  measure <- vector();
  for ( dr in dirList ) {
    directory <- sprintf("%s/Binding_sites_#1_ch_0", dr);
    fileName <- dir(directory, full.names = TRUE)[1];
    dfm <- getCoordinatesData(fileName);
    if ( funct == "Density" ) {
      calc <- match.fun(intraClusterDensity);
    }
    else if ( funct == "SSDist" ) {
      calc <- match.fun(intraClusterSS);
    }
    else if ( funct == "MPD" ) {
      calc <- match.fun(intraClusterMPD);
    }
    mes <- calc(dfm);
    measure <- append(measure, mes);
  }
  return(measure);
}


boxplotClusterMeasures <- function(dirListPer, dirListRnd, funct = c("SSDist", "Density", "MPD")[1], jitt = FALSE, ...) {
# Plots in a boxplot the intra-cluster measure that is specified by funct.
  icPer <- getIntraCusterMeasure(dirListPer, funct);
  icRnd <- getIntraCusterMeasure(dirListRnd, funct);
  lenP <- length(icPer);
  lenR <- length(icRnd);
  dfm <- data.frame(Arrangement = factor(c(rep("Periodic", lenP), rep("Random", lenR))), intraClust = c(icPer, icRnd));
  p <- ggplot(dfm, aes(Arrangement, intraClust));
  p <- p + geom_boxplot(aes(fill = Arrangement));
  if ( jitt == TRUE ) {
    p <- p + geom_jitter();
  }
  # Just to get the right legend.
  if ( funct == "Density" ) {
    p <- p + labs(y = "Intra-Cluster Density")
  }
  else if ( funct == "SSDist" ) {
    p <- p + labs(y = "Intra-Cluster Sum Squared Dist")
  }
  else if ( funct == "MPD" ) {
    p <- p + labs(y = "Intra-Cluster Mean Pair-wise Dist")
  }
  # Perform a Man-Witney-Wilcoxon test (U-test)
  uTest <- wilcox.test(icPer, icRnd);
  p <- p + ggtitle(sprintf("U-Test p value: %.3e", uTest$p.value))
  return(p)
}


clusterPairsFreqs <- function(dirList) {
# Return a list of two data frames of the contacts coocurence and its frequencies.
  # Firstly  generate a zero dfm with the appropriate names.
  dd <- getCoordinatesData(dir(sprintf("%s/Binding_sites_#1_ch_0", dirList[2]), full.names = TRUE)[1]);
  xx <- replicate(length(rownames(dd)), rep(0.0, length(rownames(dd))));
  dff <- data.frame(xx , row.names = rownames(dd));
  names(dff) <- rownames(dd);
  # The actual calculation.
  for ( dr in dirList ) {
    directory <- sprintf("%s/Binding_sites_#1_ch_0", dr);
    fileName <- dir(directory, full.names = TRUE)[1];
    dfm <- getCoordinatesData(fileName);
    clust <- dbscanClustering(dfm)$clusters;
    for ( i in 1:length(clust) ) {
      ni <- names(clust[i]);
      for ( j in i:length(clust) ) {
        nj <- names(clust[j]);
        if ( i != j && clust[i] == clust[j] ) {
          # Means that the two points belong to the same cluster so increase the frequancy by one.
          dff[ni, nj] <- dff[ni, nj] + 1 ;
          # Also populate the diametric cell of the matrix.
          dff[nj, ni] <- dff[nj, ni] + 1 ;
        }
      }
    }
  }
  # Divide with the number of experiments to get the frequency.
  ne <- length(dirList);
  freq <- dff / ne;
  # Report both.
  return(list(contacts = dff, freqs = freq))
}


clusterPairsCooccurenceMeasures <- function(dff, type = c("median", "avg", "max")[1], ....) {
# Multivalent function that return a cluster co-ocurence measure foraech interaction site.
  coOccM <- vector();
  for ( r in rownames(dff$freqs) ) {
    if ( type == "median" ) {
      mes <- median(dff$freqs[[r]])
    }
    else if ( type == "avg" ) {
      mes <- mean(dff$freqs[[r]])
    }
    else if ( type == "max" ) {
      mes <- max(dff$freqs[[r]])
    }
    coOccM <- append(coOccM, mes);
  }
  names(coOccM) <- rownames(dff$freqs);
  return(coOccM)
}


clusterPairsSpecificities <- function(dff, type = c("sum", "avg")[1], ...) {
# Return a clustering specificity index foreach interaction site. Calculated in two means one with the averages one with the sums.
  n <- length(dff$freqs);
  clSpecAVG <- vector();
  clSpecSUM <- vector();
  for ( r in rownames(dff$freqs) ) {
    coMax <- max(dff$freqs[r,]);
    cM <- which(dff$freqs[r,] == coMax)[1];
    mRow <- mean(as.numeric(dff$freqs[r,]));
    mCol <- mean(as.numeric(dff$freqs[,cM]));
    sRow <- sum(dff$freqs[r,]);
    sCol <- sum(dff$freqs[,cM]);
    csA <- (coMax) / (mCol + mRow);
    csS <- (coMax) / (sCol + sRow);
    clSpecAVG <- append(clSpecAVG, csA);
    clSpecSUM <- append(clSpecSUM, csS);
  }
  names(clSpecAVG) <- rownames(dff$freqs);
  names(clSpecSUM) <- rownames(dff$freqs);
  if ( type == "avg" ) {
    return(clSpecAVG)
  }
  else if (type == "sum" ) {
    return(clSpecSUM)
  }
}


contactFrequenciesMatrix <- function(dirList, contactRadious, type = c("circular", "linear")[1], ...) {
# Function to generate a contact frequencies matrix, equivalent (and thus easy to compare) with the matrices of 3C experiments. Calculates the contact probability of cylinders.
  # Initialise the contact probability matrix.
  dr1 <- sprintf("%s/Coordonnees_ADN_ch_0", dirList[1])
  fln <- dir(dr1, full.names = TRUE)[1];
  dd <- getCoordinatesData(fln);
  noCyl <- nrow(dd) - 1;
  # Knowing the number of cylinders let us preconstruct the contact frequency matrix.
  contaFreq <- matrix(0, nrow = noCyl, ncol = noCyl, byrow = TRUE);
  # Iterate over the experiments.
  for ( dr in dirList ) {
    directory <- sprintf("%s/Coordonnees_ADN_ch_0", dr);
    fileName <- dir(directory, full.names = TRUE)[1];
    dd <- getCoordinatesData(fileName);
    # We keep all the cylinders apart for the last one which is a duplicate in circular chromosomes.
    dd <- dd[-nrow(dd),];
    ddDist <- dist(dd, diag = TRUE, upper = TRUE);
    ddDist1 <- (as.matrix(ddDist) < contactRadious & as.matrix(ddDist) != 0) * 1; # Amazing trick!!!
    contaFreq <- contaFreq + ddDist1;
  }
  # Divide with the number of experiments to give the frequency.
  return(list(contacts = contaFreq, freqs = contaFreq / length(dirList)))
}



###############################################################################
# Higher level functions - provide plots and stats for collections of output files.
###############################################################################

# Collection fiber visualisation functions.

plotCollectionFiber <- function(directoryName, ...) {
# Plot the fibre 3D conformation iteratively.
  filesList <- dir(directoryName, full.names = TRUE);
  tl <- substr(directoryName, 1, nchar(directoryName)-1);
  tl <- gsub("_", " ", tl);
  open3d(windowRect = c(windowRect = c(0, 0, 800, 800)));
  for (fl in filesList[2:length(filesList)]) {
    plotFibre(fl, subtitle = tl, op = FALSE, ...);
    readline("Press <Enter> for the next plot.");
  }
}


plotCollectionFiberInteractions <- function(tracesDir, ...) {
# Plot the fibre 3D structure as well as the interactions iteratively.
  stop("Not implemented yet!");
}


plotCollectionTraces <- function(directoryListPer, directoryListRnd, funct, type, ...) {
# Plot the simulation traces (either PMD or Girth) from a collection of coordinates files during the whole simulation.
  if (funct == "pmd") {
    st <- "Mean Pairwise Distance"
  }
  else if (funct == "girth") {
    st <- "Gyration Radius"
  }
  else { stop(); }
  tl <- sprintf(" %s Traces", st);
  tl <- gsub("_", " ", tl);  
  perCollectedData <- getMultipleTraceData(directoryListPer, funct, type);
  rndCollectedData <- getMultipleTraceData(directoryListRnd, funct, type);
  minQuant <- min(perCollectedData[["min"]], rndCollectedData[["min"]]);
  maxQuant <- max(perCollectedData[["max"]], rndCollectedData[["max"]]); 
  ylim <- c(0.95*minQuant, maxQuant);
  xlim <- c(0, max(perCollectedData[["data"]][[1]]$mcSteps) + 1); # Check the complicated data structure list of lists.
  plot(x = xlim, y = ylim, type = "n", xlab = "MC Steps", ylab = st, main = tl, ...);
  legend("topright", legend= c("Periodic", "Random"), fill = c("blue", "red"), ...)
  for (dp in perCollectedData[["data"]]) {
    lines(dp$mcSteps, dp[[funct]], col = "blue", ...);
  }
  for (dp in rndCollectedData[["data"]]) {
    lines(dp$mcSteps, dp[[funct]], col = "red", ...);
  }
}


plotCollectionTraceRatios <- function(directoryListPer, directoryListRnd, lineW = 0.6, ...) {
# Plot the ratios between the PMD of the interaction sites and the Girth of the polymer.
  # Some elementary testing
  if (length(directoryListPer) != length(directoryListRnd)) {
    stop("The two input lists differ in size no ratio can be calculated");
  }
  perCollectedDataPMD <- getMultipleTraceData(directoryListPer, "pmd", "binding");
  rndCollectedDataPMD <- getMultipleTraceData(directoryListRnd, "pmd", "binding");
  perCollectedDataGirth <- getMultipleTraceData(directoryListPer, "girth", "coordinates");
  rndCollectedDataGirth <- getMultipleTraceData(directoryListRnd, "girth", "coordinates");
  ratiosPer <- getRatios(perCollectedDataPMD, perCollectedDataGirth);
  ratiosRnd <- getRatios(rndCollectedDataPMD, rndCollectedDataGirth);
  minRatios <- min(ratiosPer$min, ratiosRnd$min);
  maxRatios <- max(ratiosPer$max, ratiosRnd$max);
  ylim <- c(0.95*minRatios, 1.05*maxRatios);
  mcSteps <- perCollectedDataPMD[["data"]][[1]]$mcSteps;
  xlim <- c(0, max(mcSteps) + 1);
  plot(x = xlim, y = ylim, type = "n", xlab = "MC Steps", main = "Ratio traces of Interaction Sites PMD over Polymer Girth", ylab = "Inter. Sites PMD / Polymer Girth", ...);
  legend("topright", legend= c("Periodic", "Random"), fill = c("blue", "red"), ...)
  for (ratio in ratiosPer[["ratioList"]]) {
    lines(mcSteps, ratio, col = "blue", lwd = lineW, ...);
  }
  for (ratio in ratiosRnd[["ratioList"]]) {
    lines(mcSteps, ratio, col = "red", lwd = lineW, ...);
  }
}


# Cluster and interactions visualisation functions.

visualiseClusterCooccurance <- function(cooccM, n = 32, tlt = "", ...) {
# Visualise the cluster coocurance matrix.
  levelplot(as.matrix(cooccM[nrow(cooccM):1,]), aspect = "iso", scales = list(x = list(rot = 90)), col.regions = rev(heat.colors(n)), ylab = "Site Position", xlab = "Site Position", main = sprintf("Cluster Co-occurence %s", tlt), axis(2, labels = rownames(cooccM)), ...);
}


getCusterCooccurenceMeasure <- function(cooccM, funct, ...) {
# Return a vector of the required cluster co-occurence measure from a cooccurence matrix.
  measure <- vector();
  if ( funct == "avg" ) {
    return(clusterPairsCooccurenceMeasures(cooccM, type = "avg"));
  }
  else if ( funct == "max" ) {
    return(clusterPairsCooccurenceMeasures(cooccM, type = "max"));
  }
  else if ( funct == "median" ) {
    return(clusterPairsCooccurenceMeasures(cooccM, type = "median"));
  }
  else if ( funct == "spec" ) {
    return(clusterPairsSpecificities(cooccM));
  }
}


barplotsClusterCoocurances <- function(cooccMPer, cooccMRnd, tlt = "", fun = c("median", "avg", "max", "spec")[1], ...) {
# Barplot the cluster co-occurence measures.
  ccP <- getCusterCooccurenceMeasure(cooccMPer, fun);
  ccR <- getCusterCooccurenceMeasure(cooccMRnd, fun);
  posa <- names(ccP);
  posb <- names(ccR);
  lena <- length(ccP);
  lenb <- length(ccR);
  dfm <- data.frame(pos = factor(c(posa, posb), levels = c(posa, posb)), Arrangment = factor(c(rep("Periodic", lena), rep("Random", lenb))), value = c(ccP, ccR));
  # Initiate a ggplot object.
  p <- ggplot(dfm, aes(x = pos, y = value, fill = Arrangment));
  p <- p + geom_bar(stat = "identity") + facet_grid(.~Arrangment, scale = "free_x", space = "free_y");
  # tilt the x axis labels.
  p <- p + theme(axis.text.x = element_text(angle = 90));
  # Set the title
  p <- p + labs(title = sprintf("Barplots of %s cluster co-occurance %s", fun, tlt), x = "Site Position", y = "Co-occurance level");
  return(p)
}


clustergramPlot <- function(X, Y, noRepeats, x.range, y.range , COL, centresPoints) {
# Function to actually make the clustergram
  plot(0,0, col = "white", xlim = x.range, ylim = y.range, axes = F, xlab = "Number of Experiments", ylab = "PCA weighted Centres of the clusters", main = c("Clustergram of the PCA-weighted Centres of clusters", "among experimental repetitions"));
  axis(side =1, at = noRepeats);
  axis(side =2);
  abline(v = 1:noRepeats, col = "grey");
  matlines(t(X), t(Y), pch = 19, col = COL, lty = 1, lwd = 1.5); 
  xx <- ldply(centresPoints, rbind);
  points(xx$y~xx$x, pch = 19, col = "red", cex = 1.1);
}


clustergramDBscan <- function(directoryList, line.width = 0.004) {
# Manage the ploting the clustergram of the DBscan results of many repetitions of each experiment.
  fileName <- dir(sprintf("%s/Binding_sites_#1_ch_0/", directoryList[1]), full.names = TRUE)[1];
  d1 <- getCoordinatesData(fileName);
  PCA.1 <- as.matrix(d1) %*% princomp(d1)$loadings[,1]; # first principal component of our data
  noSites <- dim(d1)[1];
  if(require(colorspace)) {
    COL <- heat_hcl(noSites)[order(PCA.1)]  # line colours
  } else {
    COL <- rainbow(noSites)[order(PCA.1)] # line colors
    warning('Please consider installing the package "colorspace" for prettier colours')
  }
  line.width <- rep(line.width, noSites)
  Y <- NULL; # Y matrix
  X <- NULL; # X matrix
  Data <- NULL;
  centresPoints <- list();
  # This loop does all the job.
  for (i in 1:length(directoryList)) {
    fileName <- dir(sprintf("%sBinding_sites_#1_ch_0/", directoryList[i]), full.names = TRUE)[1];
    dd <- getCoordinatesData(fileName);
    Data <- rbind(Data, dd); # Keep the data so that you can calculate the Y.range at the end.
    clustering <- dbscanClustering(dd);
    clusterVector <- clustering$cluster;
    centres <- clustering$weighted.centres;
    # The trick to include noise is taken from the original clustergram code.
    noise <- unlist(tapply(line.width, clusterVector, cumsum))[order(seq_along(clusterVector)[order(clusterVector)])];
    y <- centres[clusterVector] + noise;
    Y <- cbind(Y, y);
    x <- rep(i, length(y));
    X <- cbind(X, x);
    centresPoints[[i]] <- data.frame(y = centres , x = rep(i, max(clusterVector)));
  }
  print(range(ldply(centresPoints, rbind)$y));
  PCA.Data <- as.matrix(Data) %*% princomp(Data)$loadings[,1];
  x.range <- range(1:length(directoryList));
  y.range <- range(floor(min(ldply(centresPoints, rbind)$y))*0.95, ceiling(max(ldply(centresPoints, rbind)$y))*1.05);
  print(y.range);
  # Plot the clustergram
  clustergramPlot(X, Y, length(directoryList), x.range, y.range, COL, centresPoints);
}


clusterCooccurenceDist <- function(contactFreqs, noCyl, molecule = c("circular", "linear")[1], includeZeros = TRUE, ...) {
# Return a dataframe of the distance (in cylinders)
## The contactFreqs matrix is a data.frame with colnames and rownames the cylinder numbers.
  dst <- vector();
  freq <- vector();
  for ( rname in rownames(contactFreqs) ) {
    for ( cname in colnames(contactFreqs) ) {
      if ( rname != cname & as.numeric(cname) > as.numeric(rname) ) { # second condition is to avoid double counts on the symmetric contactFreqs matrix.
        freq <- append(freq, contactFreqs[rname, ][[cname]]);
        if ( molecule == "linear" ) {
          dst <- append(dst, abs(as.numeric(rname) - as.numeric(cname)));
        }
        else if ( molecule == "circular" ) {
          dstMod <- min(abs(as.numeric(rname) - as.numeric(cname)), abs(noCyl - abs(as.numeric(rname) - as.numeric(cname))));
          dst <- append(dst, dstMod);
        }
      }
    }
  }
  return(data.frame(dst = dst, freq = freq))
}


plotClustCooccurDistance <- function(cooccDist, tlt = "", k = 11, type = c("contour", "smooth", "ggplot")[2], ...) {
# Plot the contact frequencies of all the pairs of interaction sites, between random and regular configurations.
  # Get some pretty colors.
  myCols <- rev(brewer.pal(k, "RdYlBu"));
  # Do te plotting with kernel contours.
  if ( type == "contour" ) {
    # Compute 2D kernel density, see MASS book, pp. 130-131
    z <- kde2d(cooccDist$dst, cooccDist$freq, n = 100);
    plot(cooccDist$dst, cooccDist$freq, pch = 19, cex = 0.33, col = "red", xlab = "Distance in cylinders", ylab = "Frequency of co-clustering", main = sprintf("Contour scatterplot for %s", tlt), ...);
    contour(z, drawlabels = FALSE, nlevels = k, col = myCols, add = TRUE);
  }
  # Do the plotting with the smoothed scatterplot.
  else if ( type == "smooth" ) {
    smoothScatter(cooccDist$dst, cooccDist$freq, nrpoints = 0.0 * length(cooccDist$dst), colramp = colorRampPalette(myCols), pch = 19, cex = .25, xlab = "Distance in cylinders", ylab = "Frequency of co-clustering", main = sprintf("Smoothed scatterplot of %s", tlt), ...);
  }
  else if ( type == "ggplot" ) {
    #TODO make a separate ggplot2 function to use facets and print both Per and Rnd.
    p <- ggplot(cooccDist, aes(dst, freq));
    p <- p + stat_bin2d(bins = ceiling(sqrt(max(cooccDist$dst)))) + labs(title = sprintf("Binned scatterplot of %s", tlt), x = "Distance in cylinders", y = "Frequency of co-clustering");
    return(p)
  }
}


visualiseContactFreqs <- function(cfm, tlt, dat = c("contacts", "freqs")[1], n = 128, square = FALSE, type = c("level", "gglot")[1], rst = FALSE, ...) {
# Generate an image plot of the contact probability matrix of the fibre.
  # cfm is the return list of the contactFrequenciesMatrix.
  if ( square == FALSE ) {
    cfmDT <- cfm$freqs;
  }
  else if ( square == TRUE ) {
    cfmDT <- apply(cfm$freqs, c(1, 2), sqrt);
  }
  levelplot(cfmDT, aspect = "iso", col.regions = rev(heat.colors(n)), useRaster = rst, ylab = "Cylinder Num.", xlab = "Cylinder Num.", main = sprintf("Contact Frequencies %s", tlt), xlim = c(1:nrow(cfm$freqs)), ylim = c(1:nrow(cfm$freqs)), ...);
}



# Extras
rotate <- function(dur = 10, ...) {
# Small short-cut to the rgl play3d
  play3d(spin3d(axis = c(1, 1, 1)), duration =  dur);
}



#############################################################
# Include some testing of the measures here... Imperative!!!!
#############################################################
