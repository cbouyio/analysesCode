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


boostrap_data <- function(expressionMatrix, bootTimes = 100, n = 10, ...) {
# Do bootstraping of the data to obtain a significance threshold.
  noGenes <- dim(expressionMatrix)[1];
  corTld  <- NULL;
  for ( i in 1:100 ) {
    randomOrder <- sample(1:noGenes, noGenes, replace = FALSE);
    newData <- expressionMatrix[randomOrder, ];
    corTot <- cor(t(newData), method = "spearman");
    corMatrix <- matrix(0, ncol = 2*n + 1, nrow = noGenes - 2*n);
    for( i in 1:(noGenes - 2*n) ) {
      corMatrix[i, ] <- corTot[i + n, i:(i + 2*n)];
    }
    corTld <- cbind(corTld, rowSums(corMatrix[, -(n + 1)]) - 1);
  }
  thresold <- quantile(as.vector(corTld), probs = c(1 - 1/500));
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

