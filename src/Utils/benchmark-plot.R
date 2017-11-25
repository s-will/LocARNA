#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)

idx <- 1
if (length(args) >= 4) {
    idx <- strtoi(args[4])
}
idx <- idx + 1

feature <- ""
if (length(args) >= 3) {
    feature <- args[3]
}

title <- ""
if (length(args) >= 2) {
    title <- args[2]
}

if (length(args) == 0) {
    cat("USAGE: benchmark-plot.R <file_name> <title> <feature_name> <index>\n")
    quit()
}

filename <- args[1]

xs <- read.table(filename)

refapsi <- gsub(".*apsi-([\\d]+).*","\\1",xs$V1,fixed=F,perl=T)
fval <- xs[[idx]]

lowess <- lowess(refapsi,fval,f=1/3)

outfilename <- paste(gsub(".tab$","",filename),"pdf",sep=".")

pdf(outfilename,width=6,height=6)
plot(refapsi,fval,xlab="APSI",ylab=feature,main=title,col=rgb(0.2,0.8,1,0.1),pch=16)
lines(lowess,lwd=3,col=rgb(0,0.2,0.8,0.8))
