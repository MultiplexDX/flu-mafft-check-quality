# R script
## © 2020, Viera Kovacova for MultiplexDX. All rights reserved. 
# sequences downloaded from gisaid.org 
# filtering in gisaid: 1.1.2018-24.6.2020, human, 4 flu types (IAV H1N1, IAV H3N2, IBV Victori, IBV Yamagata); 3 segments (PB1, PB2, PA)
# Once the influenza datasets were filtered for quality (too gappy regions or too many ambiguous bases) 
# and redundant sequences (the presence of the same clone multiple times), we checked the conservancy towards the segments.
# The decision criteria for selecting a gene segment to design a novel primer set: a selected segment must have at 
# least 200 bp long region (the length closer to 200 bp is better) with at least three conserved loci of minimal size 30 bp.

### Coputing and plotting the level of conservancy ###
# The level of the conservancy is judged using a sliding window analysis. 
# We computed the percentage of different bases than the one in the reference and got the mean of 30 consequent values. 
# Then we shifted the window by ten nucleotides. 
library(Biostrings)
library(zoo)

## IAV H1N1
# PB2
PB2_H1N1 <- readDNAMultipleAlignment("~/Documents/flu/IAV_H1N1/H1N1_PB2_segm1_filter1_mafft.fa", format = "fasta")
PSSM_PB2_H1N1 <- consensusMatrix(PB2_H1N1, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_H1N1, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_H1N1)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_H1N1[,i] == max(PSSM_PB2_H1N1[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_H1N1[,i] == max(PSSM_PB2_H1N1[,i])))), collapse = "_"))
}
     
par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_H1N1)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV H1N1 PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_H1N1_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_H1N1_slide )), y=PB2_H1N1_slide , col="cornflowerblue", pch=21, type="l",  main="IAV H1N1 PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_H1N1 <- readDNAMultipleAlignment("~/Documents/flu/IAV_H1N1/H1N1_PB1_segm2_filter1_mafft.fa", format = "fasta")
PSSM_PB1_H1N1 <- consensusMatrix(PB1_H1N1, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_H1N1, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_H1N1)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_H1N1[,i] == max(PSSM_PB1_H1N1[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_H1N1[,i] == max(PSSM_PB1_H1N1[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_H1N1)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV H1N1 PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_H1N1_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_H1N1_slide )), y=PB1_H1N1_slide , col="cornflowerblue", pch=21, type="l",  main="IAV H1N1 PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_H1N1 <- readDNAMultipleAlignment("~/Documents/flu/IAV_H1N1/H1N1_PA_segm3_filter1_mafft.fa", format = "fasta")
PSSM_PA_H1N1 <- consensusMatrix(PA_H1N1, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_H1N1, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_H1N1)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_H1N1[,i] == max(PSSM_PA_H1N1[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_H1N1[,i] == max(PSSM_PA_H1N1[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_H1N1)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV H1N1 PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_H1N1_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_H1N1_slide )), y=PA_H1N1_slide , col="cornflowerblue", pch=21, type="l",  main="IAV H1N1 PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

## IAV H3N2
# PB2
PB2_H3N2 <- readDNAMultipleAlignment("~/Documents/flu/IAV_H3N2/H3N2_PB2_segm1_filter1_mafft.fa", format = "fasta")
PSSM_PB2_H3N2 <- consensusMatrix(PB2_H3N2, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_H3N2, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_H3N2)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_H3N2[,i] == max(PSSM_PB2_H3N2[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_H3N2[,i] == max(PSSM_PB2_H3N2[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_H3N2)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV H3N2 PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_H3N2_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_H3N2_slide )), y=PB2_H3N2_slide , col="cornflowerblue", pch=21, type="l",  main="IAV H3N2 PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_H3N2 <- readDNAMultipleAlignment("~/Documents/flu/IAV_H3N2/H3N2_PB1_segm2_filter1_mafft.fa", format = "fasta")
PSSM_PB1_H3N2 <- consensusMatrix(PB1_H3N2, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_H3N2, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_H3N2)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_H3N2[,i] == max(PSSM_PB1_H3N2[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_H3N2[,i] == max(PSSM_PB1_H3N2[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_H3N2)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV H3N2 PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_H3N2_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_H3N2_slide )), y=PB1_H3N2_slide , col="cornflowerblue", pch=21, type="l",  main="IAV H3N2 PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_H3N2 <- readDNAMultipleAlignment("~/Documents/flu/IAV_H3N2/H3N2_PA_segm3_filter1_mafft.fa", format = "fasta")
PSSM_PA_H3N2 <- consensusMatrix(PA_H3N2, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_H3N2, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_H3N2)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_H3N2[,i] == max(PSSM_PA_H3N2[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_H3N2[,i] == max(PSSM_PA_H3N2[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_H3N2)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV H3N2 PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_H3N2_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_H3N2_slide )), y=PA_H3N2_slide , col="cornflowerblue", pch=21, type="l",  main="IAV H3N2 PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

## IBV Victoria
# PB2
PB2_Vict <- readDNAMultipleAlignment("~/Documents/flu/IBV_Victoria/Vict_PB2_segm1_filter1_mafft.fa", format = "fasta")
PSSM_PB2_Vict <- consensusMatrix(PB2_Vict, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_Vict, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_Vict)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_Vict[,i] == max(PSSM_PB2_Vict[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_Vict[,i] == max(PSSM_PB2_Vict[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_Vict)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV Victoria PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_Vict_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_Vict_slide )), y=PB2_Vict_slide , col="cornflowerblue", pch=21, type="l",  main="IBV Victoria PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_Vict <- readDNAMultipleAlignment("~/Documents/flu/IBV_Victoria/Vict_PB1_segm2_filter1_mafft.fa", format = "fasta")
PSSM_PB1_Vict <- consensusMatrix(PB1_Vict, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_Vict, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_Vict)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_Vict[,i] == max(PSSM_PB1_Vict[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_Vict[,i] == max(PSSM_PB1_Vict[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_Vict)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV Victoria PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_Vict_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_Vict_slide )), y=PB1_Vict_slide , col="cornflowerblue", pch=21, type="l",  main="IBV Victoria PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_Vict <- readDNAMultipleAlignment("~/Documents/flu/IBV_Victoria/Vict_PA_segm3_filter1_mafft.fa", format = "fasta")
PSSM_PA_Vict <- consensusMatrix(PA_Vict, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_Vict, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_Vict)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_Vict[,i] == max(PSSM_PA_Vict[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_Vict[,i] == max(PSSM_PA_Vict[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_Vict)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV Victoria PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_Vict_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_Vict_slide )), y=PA_Vict_slide , col="cornflowerblue", pch=21, type="l",  main="IBV Victoria PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

## IBV Yamagata
# PB2
PB2_Yama <- readDNAMultipleAlignment("~/Documents/flu/IBV_Yamagata/Yama_PB2_segm1_filter1_mafft.fa", format = "fasta")
PSSM_PB2_Yama <- consensusMatrix(PB2_Yama, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_Yama, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_Yama)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_Yama[,i] == max(PSSM_PB2_Yama[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_Yama[,i] == max(PSSM_PB2_Yama[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_Yama)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV Yamagata PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_Yama_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_Yama_slide )), y=PB2_Yama_slide , col="cornflowerblue", pch=21, type="l",  main="IBV Yamagata PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_Yama <- readDNAMultipleAlignment("~/Documents/flu/IBV_Yamagata/Yama_PB1_segm2_filter1_mafft.fa", format = "fasta")
PSSM_PB1_Yama <- consensusMatrix(PB1_Yama, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_Yama, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_Yama)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_Yama[,i] == max(PSSM_PB1_Yama[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_Yama[,i] == max(PSSM_PB1_Yama[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_Yama)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV Yamagata PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_Yama_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_Yama_slide )), y=PB1_Yama_slide , col="cornflowerblue", pch=21, type="l",  main="IBV Yamagata PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_Yama <- readDNAMultipleAlignment("~/Documents/flu/IBV_Yamagata/Yama_PA_segm3_filter1_mafft.fa", format = "fasta")
PSSM_PA_Yama <- consensusMatrix(PA_Yama, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_Yama, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_Yama)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_Yama[,i] == max(PSSM_PA_Yama[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_Yama[,i] == max(PSSM_PA_Yama[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_Yama)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV Yamagata PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_Yama_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_Yama_slide )), y=PA_Yama_slide , col="cornflowerblue", pch=21, type="l",  main="IBV Yamagata PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")


## IAV combo
# PB2
PB2_IAV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV/Cons95_PB2_segm1_H1N1_H3N2_mafft.fa", format = "fasta")
PSSM_PB2_IAV <- consensusMatrix(PB2_IAV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_IAV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_IAV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_IAV[,i] == max(PSSM_PB2_IAV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_IAV[,i] == max(PSSM_PB2_IAV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_IAV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV consensus PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_IAV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_IAV_slide )), y=PB2_IAV_slide , col="cornflowerblue", pch=21, type="l",  main="IAV consensus PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_IAV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV/Cons95_PB1_segm2_H1N1_H3N2_mafft.fa", format = "fasta")
PSSM_PB1_IAV <- consensusMatrix(PB1_IAV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_IAV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_IAV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_IAV[,i] == max(PSSM_PB1_IAV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_IAV[,i] == max(PSSM_PB1_IAV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_IAV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV consensus PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_IAV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_IAV_slide )), y=PB1_IAV_slide , col="cornflowerblue", pch=21, type="l",  main="IAV consensus PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_IAV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV/Cons95_PA_segm3_H1N1_H3N2_mafft.fa", format = "fasta")
PSSM_PA_IAV <- consensusMatrix(PA_IAV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_IAV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_IAV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_IAV[,i] == max(PSSM_PA_IAV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_IAV[,i] == max(PSSM_PA_IAV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_IAV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV consensus PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_IAV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_IAV_slide )), y=PA_IAV_slide , col="cornflowerblue", pch=21, type="l",  main="IAV consensus PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

## IBV combo
# PB2
PB2_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IBV/Cons95_PB2_segm1_Victoria_Yamagata_mafft.fa", format = "fasta")
PSSM_PB2_IBV <- consensusMatrix(PB2_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_IBV[,i] == max(PSSM_PB2_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_IBV[,i] == max(PSSM_PB2_IBV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_IBV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV consensus PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_IBV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_IBV_slide )), y=PB2_IBV_slide , col="cornflowerblue", pch=21, type="l",  main="IBV consensus PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IBV/Cons95_PB1_segm2_Victoria_Yamagata_mafft.fa", format = "fasta")
PSSM_PB1_IBV <- consensusMatrix(PB1_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_IBV[,i] == max(PSSM_PB1_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_IBV[,i] == max(PSSM_PB1_IBV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_IBV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV consensus PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_IBV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_IBV_slide )), y=PB1_IBV_slide , col="cornflowerblue", pch=21, type="l",  main="IBV consensus PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IBV/Cons95_PA_segm3_Victoria_Yamagata_mafft.fa", format = "fasta")
PSSM_PA_IBV <- consensusMatrix(PA_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_IBV[,i] == max(PSSM_PA_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_IBV[,i] == max(PSSM_PA_IBV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_IBV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IBV consensus PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_IBV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_IBV_slide )), y=PA_IBV_slide , col="cornflowerblue", pch=21, type="l",  main="IBV consensus PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# IAV and IVB combo
# PB2
PB2_IAV_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV_IBV/Cons95_PB2_segm1_H1N1_H3N2_Vict_Yama_mafft.fa", format = "fasta")
PSSM_PB2_IAV_IBV <- consensusMatrix(PB2_IAV_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB2_IAV_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB2_IAV_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB2_IAV_IBV[,i] == max(PSSM_PB2_IAV_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB2_IAV_IBV[,i] == max(PSSM_PB2_IAV_IBV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB2_IAV_IBV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV and IBV consensus PB2 Segment1\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB2_IAV_IBV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB2_IAV_IBV_slide )), y=PB2_IAV_IBV_slide , col="cornflowerblue", pch=21, type="l",  main="IAV and IBV consensus PB2 Segment1\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PB1
PB1_IAV_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV_IBV/Cons95_PB1_segm2_H1N1_H3N2_Vict_Yama_mafft.fa", format = "fasta")
PSSM_PB1_IAV_IBV <- consensusMatrix(PB1_IAV_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_IAV_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_IAV_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_IAV_IBV[,i] == max(PSSM_PB1_IAV_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_IAV_IBV[,i] == max(PSSM_PB1_IAV_IBV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PB1_IAV_IBV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV and IBV consensus PB1 Segment2\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PB1_IAV_IBV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PB1_IAV_IBV_slide )), y=PB1_IAV_IBV_slide , col="cornflowerblue", pch=21, type="l",  main="IAV and IBV consensus PB1 Segment2\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

# PA
PA_IAV_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV_IBV/Cons95_PA_segm3_H1N1_H3N2_Vict_Yama_mafft.fa", format = "fasta")
PSSM_PA_IAV_IBV <- consensusMatrix(PA_IAV_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_IAV_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_IAV_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_IAV_IBV[,i] == max(PSSM_PA_IAV_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_IAV_IBV[,i] == max(PSSM_PA_IAV_IBV[,i])))), collapse = "_"))
}

par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=0.85, cex.lab=1)
plot(x=c(1:dim(PSSM_PA_IAV_IBV)[2]), y=maxProb, col="cornflowerblue", pch=".", type="p", main="IAV and IBV consensus PA Segment3\nThe highest probability of the assigned base per position", ylab="max probability", xlab="position (bp)")
PA_IAV_IBV_slide <- rollapply(maxProb, width=30, by=10, mean)
plot(x=c(1:length(PA_IAV_IBV_slide )), y=PA_IAV_IBV_slide , col="cornflowerblue", pch=21, type="l",  main="IAV and IBV consensus PA Segment3\nThe highest probability of the assigned base per position", ylab="averaged max probability", xlab="sliding window (width=30bp, shift by 10 bp)")

### seqLogo
library(Biostrings)
library(seqLogo)

# IAV
PB1_IAV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IAV/Cons95_PB1_segm2_H1N1_H3N2_mafft.fa", format = "fasta")
PSSM_PB1_IAV <- consensusMatrix(PB1_IAV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PB1_IAV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PB1_IAV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PB1_IAV[,i] == max(PSSM_PB1_IAV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PB1_IAV[,i] == max(PSSM_PB1_IAV[,i])))), collapse = "_"))
}
# PB1, segment2 - set2
# 1900-2205: 10-21-237-21-61-22
# 2157-2288
par(mar=c(3,3,2.5,0.5), mgp=c(2, 0.75,0), cex.main=0.9, cex.axis=0.85, cex.lab=1)
# 1890:1940,2145:2205
plot(x=c(1:111), y=maxProb[c(1891:1930,2147:2217)], type="p", pch=19, col="cornflowerblue", xaxt="na", xlab="", ylab="max probability", main="IAV consensus - PB1, segment2 - RT-qPCR set 2 (product length 306 bp)")
rect(10, 0.5, 29, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 19, y=0.8, "forward\nprimer", cex=1.5, srt=0, col="gray50")
axis(1, at = c(1:111) , labels = myChars[c(1891:1930,2147:2217)], cex.axis=0.5, las=2)
abline(v=40.25, col="coral2", lty=2)
abline(v=40.75, col="coral2", lty=2)
text(x = 39.5, y=0.7, "216 skipped bases", cex=0.75, srt=90)
rect(52, 0.5, 74, 1, col = "cornsilk3", border = FALSE, density = 25, angle = 45)
text(x = 62, y=0.8, "probe", cex=1.5, srt=0, col="gray50")
rect(82, 0.5, 99, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 90, y=0.8, "reverse\nprimer", cex=1.5, srt=0, col="gray50")
# IAVconsensus_set2_probabilityPlot

# PB1, segment2 - set1
# 2157-2288: 
par(mar=c(3,3,2.5,0.5), mgp=c(2, 0.75,0), cex.main=0.9, cex.axis=0.85, cex.lab=1)
# 2147:2215, 2257:2298
plot(x=c(1:111), y=maxProb[c(2147:2215, 2257:2298)], type="p", pch=19, col="cornflowerblue", xaxt="na", xlab="", ylab="max probability", main="IAV consensus - PB1, segment2 - RT-qPCR set 1 (product length 132 bp)")
rect(11, 0.5, 31, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 19, y=0.8, "forward\nprimer", cex=1.5, srt=0, col="gray50")
axis(1, at = c(1:111) , labels = myChars[c(2147:2215, 2257:2298)], cex.axis=0.5, las=2)
rect(39, 0.5, 59, 1, col = "cornsilk3", border = FALSE, density = 25, angle = 45)
text(x = 46, y=0.8, "probe", cex=1.5, srt=0, col="gray50")
text(x = 68.5, y=0.7, "41 skipped bases", cex=0.75, srt=90)
abline(v=69.25, col="coral2", lty=2)
abline(v=69.75, col="coral2", lty=2)
rect(80, 0.5, 101, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 90, y=0.8, "reverse\nprimer", cex=1.5, srt=0, col="gray50")
# IAVconsensus_set1_probabilityPlot

# IBV
PA_IBV <- readDNAMultipleAlignment("~/Documents/flu/IAV_IBV_consensus/IBV/Cons95_PA_segm3_Victoria_Yamagata_mafft.fa", format = "fasta")
PSSM_PA_IBV <- consensusMatrix(PA_IBV, as.prob = TRUE, baseOnly=TRUE)
maxProb <- apply(PSSM_PA_IBV, 2, max)
myCols <- c()
myChars <- c()
for (i in c(1:dim(PSSM_PA_IBV)[2])){
  myCols <- c(myCols, unlist(which(PSSM_PA_IBV[,i] == max(PSSM_PA_IBV[,i])))[1])
  myChars <- c(myChars, paste0(unlist(names(which(PSSM_PA_IBV[,i] == max(PSSM_PA_IBV[,i])))), collapse = "_"))
}
# PA, segment3 - set1
# 1423-1669
par(mar=c(3,3,2.5,0.5), mgp=c(2, 0.75,0), cex.main=0.9, cex.axis=0.85, cex.lab=1)
plot(x=c(1:112), y=maxProb[c(1423:1462,1598:1669)], type="p", pch=19, col="cornflowerblue", xaxt="na", xlab="", ylab="max probability", main="IBV consensus - PA, segment3 - RT-qPCR set 1 (product length 227 bp)")
rect(11, 0.5, 30, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 18, y=0.8, "forward\nprimer", cex=1.5, srt=0, col="gray50")
axis(1, at = c(1:112) , labels = myChars[c(1423:1462,1598:1669)], cex.axis=0.5, las=2)
abline(v=40.25, col="coral2", lty=2)
abline(v=40.75, col="coral2", lty=2)
text(x = 39.5, y=0.7, "135 skipped bases", cex=0.75, srt=90)
rect(51, 0.5, 73, 1, col = "cornsilk3", border = FALSE, density = 25, angle = 45)
text(x = 61, y=0.8, "probe", cex=1.5, srt=0, col="gray50")
rect(81, 0.5, 102, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 90, y=0.8, "reverse\nprimer", cex=1.5, srt=0, col="gray50")
# IBVconsensus_set1_probabilityPlot

# PA, segment3 - set2
# 1600-1793
length(c(1600:1639,1713:1793))
par(mar=c(3,3,2.5,0.5), mgp=c(2, 0.75,0), cex.main=0.9, cex.axis=0.85, cex.lab=1)
plot(x=c(1:121), y=maxProb[c(1600:1639,1713:1793)], type="p", pch=19, col="cornflowerblue", xaxt="na", xlab="", ylab="max probability", main="IBV consensus - PA, segment3 - RT-qPCR set 2 (product length 174 bp)")
rect(11, 0.5, 30, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 18, y=0.8, "forward\nprimer", cex=1.5, srt=0, col="gray50")
axis(1, at = c(1:121) , labels = myChars[c(1600:1639,1713:1793)], cex.axis=0.5, las=2)
abline(v=40.25, col="coral2", lty=2)
abline(v=40.75, col="coral2", lty=2)
text(x = 39.5, y=0.7, "73 skipped bases", cex=0.75, srt=90)
rect(51, 0.5, 76, 1, col = "cornsilk3", border = FALSE, density = 25, angle = 45)
text(x = 62, y=0.8, "probe", cex=1.5, srt=0, col="gray50")
rect(88, 0.5, 111, 1, col = "bisque", border = FALSE, density = 25, angle = 45)
text(x = 99, y=0.8, "reverse\nprimer", cex=1.5, srt=0, col="gray50")
# IBVconsensus_set2_probabilityPlot

