# 3D plotting and statistics for the folding program output.

library("rgl");
library("stringr");
library("fpc");
library("plyr");
library("MASS");
library("RColorBrewer");
library("ggplot2");


########################
## Low level functions - Access, store and basic calculations on folding data.
########################

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


intraClusterSS <- function(coordsFile, e = 90, cp = 2, ...) {
# Return a vector (of length equal to the number of clusters) of the intra-cluster sum of squares of each cluster.
  dframe <- getCoordinatesData(coordsFile);
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


intraClusterDensity <- function(coordsFile, e = 90, cp = 2, ...) {
# Return a vector (of length equal to the number of clusters) of the density of each cluster.
  dframe <- getCoordinatesData(coordsFile);
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


intraClusterMPD <- function(coordsFile, e = 90, cp = 2, ...) {
# Return a vector (of length equal to the number of clusters) of the MPD of each cluster.
  dframe <- getCoordinatesData(coordsFile);
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


jaccardIndex <- function(labelsA, labelsB, ...) {
# Return the Jaccard similarity index between two sets.
  u <- union(labelsA, labelsB);
  i <- intersect(labelsA, labelsB);
  return(i/u) # The Jaccard similarity index
}


clusterPairsFreqs <- function(dirList) {
# Return a matrix (data frame) of the contacts coocurence.
  # Firstly  generate a zero dfm with the appropriate names.
  dd <- getCoordinatesData(dir(sprintf("%s/Binding_sites_#1_ch_0", dirList[1]), full.names = TRUE)[1]);
  xx <- replicate(length(rownames(dd)), rep(0, length(rownames(dd))));
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
          dff[ni, ][[nj]] <- dff[ni, ][[nj]] + 1;
        }
      }
    }
  }
  # Divide with the number of experiments to get the frequency.
  ne <- length(dirList);
  return(dff / ne)
}


clusterPairsCooccurence <- function(dff, ...) {
# Return the max of cluster co-occurance foreach interaction site.
  clustContact <- vector();
  for ( r in rownames(dff) ) {
    coocc <- max(dff[r,], t(dff)[,r]);
    clustContact <- append(clustContact, coocc);
  }
  names(clustContact) <- rownames(dff);
  return(clustContact)
}


clusterSiteSpecificity <- function(dff, ...) {
# Return the average of easch line of the matrix dff.
#  n <- length(dff);
#  for ( r in )
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
    directory <- sprintf("%s%s", dr, suff);
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



#######################
## High level functions - Analyse and plot chromosome folding single file data.
#######################

plotFiber <- function(fileName, lineWidth = 2.5, subtitle = "Demo", op = TRUE, ...) {
# Baseline 3d plotting function of the fibre.
  if (op == TRUE) {
    open3d(windowRect = c(windowRect = c(0, 0, 800, 800)));
  }
  df <- getCoordinatesData(fileName);
  mcStep <- as.numeric(str_extract(fileName, "[0-9]{11}"));
  plot3d(df, type = "l", lwd = lineWidth, xlab = "x or 1", ylab = "y or 2", zlab = "z or 3", box = FALSE, top = TRUE, main = sprintf("MC Step: %s", mcStep), sub = subtitle, cex = 0.75, col = "darkgray", ...);
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
    text3d(dfInter$rfin.0. + 10, dfInter$rfin.1. + 10, dfInter$rfin.2. + 10, texts = rownames(dfInter), ...);
  }
}


makeBoxplots <- function(outDirListPer, outDirListRnd, funct, type, colours = c("turquoise", "tomato"), ...) {
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
  boxplot(perData, rndData, varwidth = TRUE, notch = TRUE, names = c("Periodic", "Random"), ylab = "Nanometres", main = sprintf("%s Boxplots", label), boxwex = 0.6, col = colours, ...);
  legend("topright", legend = sprintf("U-test statistic: %.2f\nU-test p value: %.4f", wilcox.test(perData, rndData)$statistic, wilcox.test(perData, rndData)$p.value))
}


plotTraces <- function(directoryName, funct, ...) {
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


clustergramDBscan <- function(directoryList, line.width = 0.004) {
# Plot the clustergram of the DBscan results of many repetitions of each experiment.
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


clustergramPlot <- function(X, Y, noRepeats, x.range, y.range , COL, centresPoints) {
# Plot the clustergram
  plot(0,0, col = "white", xlim = x.range, ylim = y.range, axes = F, xlab = "Number of Experiments", ylab = "PCA weighted Centres of the clusters", main = c("Clustergram of the PCA-weighted Centres of clusters", "among experimental repetitions"));
  axis(side =1, at = noRepeats);
  axis(side =2);
  abline(v = 1:noRepeats, col = "grey");
 
  matlines(t(X), t(Y), pch = 19, col = COL, lty = 1, lwd = 1.5);
 
  xx <- ldply(centresPoints, rbind);
  points(xx$y~xx$x, pch = 19, col = "red", cex = 1.1);
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


plotContactFrequencies <- function() {
#
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
