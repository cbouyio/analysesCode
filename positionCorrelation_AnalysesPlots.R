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
  return(sort(modCoords))
}


neighbour_correlation <- function(expressionMatrix, genePositions, n = 10, cc = c("spearman", "pearson", "kendall")[1], method = c("1D" ,"3D")[1], per = FALSE, geneList = FALSE, cyclic = TRUE, ...) {
# Calculate the correlation coefficient between the gene expression profiles of a gene and the average of its N neighbours.
  # Read files and generate the data structures.
  geMat <- read.table(expressionMatrix, header = TRUE);
  gp <- read.table(genePositions);
  gpp <- gp[[2]];
  names(gpp) <- gp[[1]];
  # Keep only the genes for which we have expression profiles from the gene expression matrix.
  genePos <- gpp[rownames(ge)];
  noGenes <- length(genePos);  
  # Transform the position to circular coordinates (if is asked).
  if ( method == "3D" ) {
    geneCoords <- modulo_coordinates(genePos, per);
  }
  else if ( method == "1D" ) {
    geneCoords <- sort(genePos);
  }
  # Calculating correlations.
  for ( i in 1:nrow(geMat) ) {
    geneExPr <- geMat[i,];
    genename <- rownames(geneExPr);
    geneIndex <- which(geneCoords == geneCoords[[geneName]]);
    # Check the positioning og the gene index and generate the neighbourhood.
    if ( geneIndex + n <= noGenes & geneIndex - n >= 1 ) {
      neighCoords <- c(seq(geneIndex - n , n), seq(geneIndex + 1, n));
    }
    else if () {
    }
    else if () {
    }
  }

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


# test
#data <- read.table("/Users/samanhasan/Desktop/gpr/1gpr(work)/TExpressionPL.txt")
Pos<- read.table("/Users/samanhasan/Desktop/gpr/1gpr(work)/TPOsitionPL.txt", header=T)

# Modulus
per <- 300295
Modulus  <-  Pos[,3] - (Pos[,3]%/%per * per);
res <-  cbind(Pos,Modulus);

### Score I: Spearman correlation value
# to obtain threshold
n=10
m=14

# correlation of real data
cor.tot <- cor(t(data), method="spearman")
cor.matrix = matrix(0, ncol=2*n+1, nrow=nbgene-2*n)
for(i in 1:(nbgene-2*n) ) cor.matrix[i,] <- cor.tot[i+n, i:(i+2*n)]
correlation =  rowSums(cor.matrix)-1

# test by 1-D
par(mfrow=c(1,1))
for(i in 1:m) {  
chri = subset(res, res[,2]==paste("chr",i,sep=""))
print(unique(chri[,2]))
pos = chri[,3]
plot(pos[(n+1):(nbgene-n)], correlation,main=paste("Chromosome",i,sep=""), ylim=c(-20,20),xlab="Genes position", ylab="TC Score", type="p",col="red")
abline(h=thresold, col="red")
}

# test by 3-D
par(mfrow=c(1,1))
for(i in 1:m) {  
  chri = subset(res, res[,2]==paste("chr",i,sep=""))
  print(unique(chri[,2]))
  pos = chri[,4]
  plot(pos[(n+1):(nbgene-n)], correlation,main=paste("Chromosome",i,sep=""), ylim=c(-20,20),xlab="Genes position", ylab="TC Score", type="p",col="red")
  abline(h=thresold, col="red")
}

