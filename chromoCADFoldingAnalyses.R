# 3D plotting and statistics for the folding program output.

library("rgl");
library("stringr");
library("fpc");
library("plyr");
library("MASS");
library("RColorBrewer");
library("ggplot2");
library("gridExtra");


########################
## Low level functions - Access, store and basic calculations on folding data.
########################

# Read data from singular files as well as collections (directories etc.).

getCoordinatesData <- function(fileName) {
# Return data as a frame with the 3D coordinates of the end-points of each cylinder.
# The cylinder number becomes the row name in the DF, and will be used to identify the point.
  d <- read.table(fileName, header = TRUE, comment.char = "&");
  rn <- as.character(d$X.mono);
  # Correction for cyclical molecules for the last cylinder which has the same coordinates with the first.
  if (rn[1] == rn[length(rn)]) {
    rn[length(rn)] <- paste(rn[length(rn)], "-c", sep = "");
  }
  coordsDF <- data.frame(d[2:4], row.names = rn);
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
    icSS <- withinClusterSS(dfCluster);
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




#######################
## High level functions - Analyse and plot chromosome folding single file data.
#######################

plotFiber <- function(fileName, lineWidth = 3, subtitle = "Demo", op = TRUE, ...) {
# Baseline 3d plotting function of the fibre.
  if (op == TRUE) {
    open3d(windowRect = c(windowRect = c(0, 0, 800, 800)));
  }
  df <- getCoordinatesData(fileName);
  mcStep <- as.numeric(str_extract(fileName, "[0-9]{11}"));
  plot3d(df, type = "l", lwd = lineWidth, xlab = "x or 1", ylab = "y or 2", zlab = "z or 3", box = FALSE, top = TRUE, col = "darkgray", ...);
  title3d(main = sprintf("MC Step: %s", mcStep), sub = subtitle);
}


plotFiberInteractions <- function(directoryName, sizeLine = 2, sizePoint = 5, ...) {
# Plot the fibre as well as the interaction sites.
# Works only for interaction site families on one chromosome!
  tl <- substr(directoryName, 1, nchar(directoryName)-1);
  tl <- gsub("_", " ", tl);
  bindingSitesDirsList <- Sys.glob(sprintf("%s/Binding_sites_#*", directoryName));
  fileNameFiber <- dir(sprintf("%s/Coordonnees_ADN_ch_0/", directoryName), full.names = TRUE);
  plotFiber(fileNameFiber, sizeLine, subtitle = "", ...);
  # Get the number of interacting families to assign different colour.
  colours <- rainbow(length(bindingSitesDirsList));
  for (i in 1:length(bindingSitesDirsList)) {
    fileNameInter <- dir(sprintf("%s/Binding_sites_#%i_ch_0/", directoryName, i), full.names = TRUE);
    dfInter <- getCoordinatesData(fileNameInter);
    points3d(dfInter$rfin.0., dfInter$rfin.1., dfInter$rfin.2., size = sizePoint, col = colours[i], ...);
    texts3d(dfInter$rfin.0., dfInter$rfin.1., dfInter$rfin.2., texts = rownames(dfInter), adj = c(-0.1,-0.1), ...);
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


plotClusterNumbers <- function(dirListPer, dirListRnd, type = c("bar", "box")[1], ...) {
# Plot the barplot of cluster numbers plus error bars.
  ncPer <- clusterNumber(dirListPer);
  ncRnd <- clusterNumber(dirListRnd);
  lenP <- length(ncPer);
  lenR <- length(ncRnd);
  dfm <- data.frame(Arrangement = factor(c(rep("Periodic", lenP), rep("Random", lenR))), NoClusters = c(ncPer, ncRnd));
  p1 <- ggplot(dfm, aes(fill = Arrangement, factor(NoClusters)));
  p1 <- p1 + geom_bar(position = "dodge") + labs(x = "Number of Clusters", y = "Count") ;
  p2 <- ggplot(dfm, aes(Arrangement, NoClusters));
  p2 <- p2 + geom_boxplot(notch = TRUE, aes(fill = Arrangement)) + geom_jitter() + labs(y = "Number of Clusters");
  if ( type == "bar" ) {
    return(p1);
  }
  else if ( type == "box" ) {
    return(p2);
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


boxplotClusterMeasures <- function(dirListPer, dirListRnd, funct = c("SSDist", "Density", "MPD")[1], ...) {
# Plots in a boxplot the intra-cluster measure that is specified by funct.
  icPer <- getIntraCusterMeasure(dirListPer, funct);
  icRnd <- getIntraCusterMeasure(dirListRnd, funct);
  lenP <- length(icPer);
  lenR <- length(icRnd);
  dfm <- data.frame(Arrangement = factor(c(rep("Periodic", lenP), rep("Random", lenR))), intraClust = c(icPer, icRnd));
  p <- ggplot(dfm, aes(Arrangement, intraClust));
  p <- p + geom_boxplot(notch = TRUE, aes(fill = Arrangement)) + geom_jitter();
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
# Return a matrix (data frame) of the contacts coocurence.
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


clusterPairsCooccurencAVG <- function(dff, ....) {
# Return the average cluster co-ocurence frequency foraech interaction site.
  coOccAvg <- vector();
  for ( r in rownames(dff$freqs) ) {
    avg <- mean(dff$freqs[r,]);
    coOccAvg <- append(coOccAvg, avg);
  }
  names(coOccAvg) <- rownames(dff$freqs);
  return(coOccAvg)
}


clusterPairsCooccurenceMax <- function(dff, ...) {
# Return the maximum cluster co-occurance frequency foreach interaction site.
  coOccMax <- vector();
  for ( r in rownames(dff$freqs) ) {
    co <- max(dff$freqs[r,]);
    coOccMax <- append(coOccMax, co);
  }
  names(coOccMax) <- rownames(dff$freqs);
  return(coOccMax)
}


clusterSiteSpecificities <- function(dff, ...) {
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
    csA <- (2 * coMax) / (mCol + mRow);
    csS <- (coMax ** 2) / (sCol + sRow);
    clSpecAVG <- append(clSpecAVG, csA);
    clSpecSUM <- append(clSpecSUM, csS);
  }
  names(clSpecAVG) <- rownames(dff$freqs);
  names(clSpecSUM) <- rownames(dff$freqs);
  return(list(avg=clSpecAVG, sum=clSpecSUM))
}


contactFreqs <- function(dirList, thresDist = 90, ...) { # Threshold distance in nanometers.
# Return a contact frequency matrix.
  # Firstly generate a zero dfm with the appropriate names.
  dd <- getCoordinatesData(dir(sprintf("%s/Coordonnees_ADN_ch_0", dirList[1]), full.names = TRUE)[1]);
  xx <- replicate(length(rownames(dd)), rep(0, length(rownames(dd))));
  cfm <- data.frame(xx , row.names = rownames(dd));
  names(cfm) <- rownames(dd);
  # Calculate the contact frequencies for cylinders that are under the defined distance (thresDist).
  for ( dr in dirList ) {
    directory <- sprintf("%s/Coordonnees_ADN_ch_0", dr);
    fileName <- dir(directory, full.names = TRUE)[1];
    dfm <- getCoordinatesData(fileName);
    for ( r in 1:nrow(dfm) ) {
      rw <- dfm[r,];
    }
  }
  ne = length(dirList);  
  return(cfm / ne)
}



########################
# Higher level functions - provide plots and stats for collections of output files.
########################

plotCollectionFiber <- function(directoryName, ...) {
# Plot the fibre 3D conformation iteratively.
  filesList <- dir(directoryName, full.names = TRUE);
  tl <- substr(directoryName, 1, nchar(directoryName)-1);
  tl <- gsub("_", " ", tl);
  open3d(windowRect = c(windowRect = c(0, 0, 800, 800)));
  for (fl in filesList[2:length(filesList)]) {
    plotFiber(fl, subtitle = tl, op = FALSE, ...);
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


plotClusterCoocurances <- function(contactFreqs, chromLen, molecule = c("circular", "linear")[1], plotType = c("contour", "smooth")[1], ...) {
# Plot the contact frequencies of all the pairs of interaction sites, between random and regular configurations.
  dst <- vector();
  freq <- vector();
  for ( rname in rownames(contactFreqs) ) {
    for ( cname in colnames(contactFreqs) ) {
      if ( rname != cname & contactFreqs[rname, ][[cname]] != 0 ) {
        freq <- append(freq, contactFreqs[rname, ][[cname]]);
        if ( molecule == "linear" ) {
          dst <- append(dst, abs(as.numeric(rname) - as.numeric(cname)));
        }
        else if ( molecule == "circular" ) {
          dstMod <- min(abs(as.numeric(rname) - as.numeric(cname)), abs(abs(chromLen - as.numeric(cname)) + as.numeric(rname)));
          dst <- append(dst, dstMod);
        }
      }
    }
  }
  # Get some pretty colors.
  k <- 11;
  myCols <- rev(brewer.pal(k, "RdYlBu"));
  # Do te plotting with kernel contours.
  if ( plotType == "contour" ) {
    # Compute 2D kernel density, see MASS book, pp. 130-131
    z <- kde2d(dst, freq, n = 100);
    plot(dst, freq, pch = 19, cex = 0.33, col = "red", xlab = "Distance in Cyliders", ylab = "Frequency of Common Clustering");
    contour(z, drawlabels = FALSE, nlevels = k, col = myCols, add = TRUE);
  }
  # Do the plotting with the smoothed scatterplot.
  else if ( plotType == "smooth" ) {
    smoothScatter(dst, freq, nrpoints = 0.0 * length(dst), colramp = colorRampPalette(myCols), pch = 19, cex = .25, xlab = "Distance in Cyliders", ylab = "Frequency of Common Clustering");
  }
}


contactFrequenciesMatrix <- function(outDirs, contactRadious, ...) {
# Function to generate a contact frequencies matrix, equivalent (and thus easy to compare) with the matrices of 3C experiments.
  return(0)
}

# Extras
rotate <- function(dur = 10, ...) {
# Small short-cut to the rgl play3d
  play3d(spin3d(axis = c(1, 1, 1)), duration =  dur);
}



#############################################################
# Include some testing of the measures here... Imperative!!!!
#############################################################
