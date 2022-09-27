#!/usr/bin/env Rscript

## benchmark plot with experimental support for several data series
##
## series must be separated by ','; the separator symbol is not allowed in input
## file names

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

outformat <- "pdf"
if (length(args) >= 5) {
    outformat <- args[5]
}

if (length(args) == 0) {
    cat("USAGE: benchmark-multi-plot.R <file_name> <title> <feature_name> <index> <out-format>\n")
    quit()
}

filenames <- unlist(strsplit(args[1],","))

outfilename <- paste(paste(title,feature,sep="-"),outformat,sep=".")

if (outformat == "pdf") {
    pdf(outfilename,width=6,height=6)
} else if (outformat == "svg") {
    svg(outfilename,width=6,height=6)
} else if (outformat == "png") {
    png(outfilename,width=6,height=6)
}

ymax = 1
if(feature=="SCI") {ymax=1.5}

plot(NULL,xlim=c(15,95),ylim=c(0,ymax),xlab="APSI",ylab=feature,main=title)

i=0
colors = c(rgb(0.2,0.8,1,0.03),rgb(1,0.2,0.4,0.03),rgb(0.4,1,0.2,0.03))
strong_colors = c(rgb(0.2,0.8,1,0.9),rgb(1,0.2,0.4,0.9),rgb(0.4,1,0.2,0.9))

for (filename in filenames) {
    i = i+1
    xs <- read.table(filename)

    refapsi <- gsub(".*apsi-([\\d]+).*","\\1",xs$V1,fixed=F,perl=T)
    fval <- xs[[idx]]

    lowess <- lowess(refapsi,fval,f=1/3)

    points(refapsi,fval,col=colors[i],pch=16)
    lines(lowess,lwd=6,col=strong_colors[i])
}
