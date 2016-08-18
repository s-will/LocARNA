gargs <- commandArgs()
t1<-read.table(gargs[5])
t2<-read.table(gargs[6])

correlations=numeric(length(t1))

for (i in 1:length(t1)) { 
   x1 <- t1[[i]];
   x2 <- t2[[i]];
   x1[i] <- 1-sum(x1);
   x2[i] <- 1-sum(x2);
   correlations[i]<-abs(cor(x1,x2));
}

print(correlations)
print(mean(correlations))

