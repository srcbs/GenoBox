# parse velvet output files and report best assembly

# input are <infile> <ksize> [<infile> <ksize> ... ]
cat("-- reading arguments\n", sep = "")
Args <- commandArgs();


# calc N50
calcN50.f = function(lengths) {
   con.v = rev(sort(lengths))
   half = sum(con.v) / 2
   csum = cumsum(con.v)
   index = sum(csum <= half)
   con.v[index]
}

parse.f = function(stats, ksize) {
   # read in stats files
   data <- read.csv(stats, header=TRUE, sep="\t")
   
   contig_lengths = data$lgth+ksize-1
   contig_sum = sum(contig_lengths)
   contig_max = max(contig_lengths)
   
   # N50
   N50 = calcN50.f(contig_lengths)
   
   # no. nodes
   nodes = nrow(data)
   
   return (c(contig_sum, N50, nodes, contig_max))
}


#Args = c("", "", "test_33/stats.txt", 33, "test_37/stats.txt", 37, "test_41/stats.txt", 41, "test_45/stats.txt", 45, "test_49/stats.txt", 49)
#Args = c("/tools/bin/R-2.12", "--vanilla", "Kleb_33/stats.txt", "33", "Kleb_37/stats.txt", "37", "Kleb_41/stats.txt", "41", "Kleb_45/stats.txt", "45", "Kleb_49/stats.txt", "49", "Kleb_53/stats.txt", "53", "Kleb_57/stats.txt", "57", "Kleb_61/stats.txt", "61", "Kleb_65/stats.txt", "65", "Kleb_69/stats.txt", "69", "Kleb_73/stats.txt", "73")
#Args = c("", "", "assembly_29/stats.txt", 29, "assembly_33/stats.txt", 33, "assembly_37/stats.txt", 37, "assembly_41/stats.txt", 41, "assembly_45/stats.txt", 45, "assembly_49/stats.txt", 49, "assembly_53/stats.txt", 53, "assembly_57/stats.txt", 57, "assembly_61/stats.txt", 61, "assembly_65/stats.txt", 65, "assembly_69/stats.txt", 69, "assembly_73/stats.txt", 73)

# parse inputs to matrix
row = 0
df = as.data.frame(matrix(NA, ncol=6, nrow=(length(Args)-2)/2))
not_finished = c()
colnames(df) = c("Ksize", "Sum", "N50", "Nodes", "MaxContig", "1/Nodes")
for (a in seq(3,length(Args),2)) {
   if (file.exists(Args[a])) {
      row = row + 1
      df[row,1] = Args[a+1]
      df[row,c(2,3,4,5)] = parse.f(Args[a], as.integer(Args[a+1]))
   } else {
      write(paste(c("Warning: assembly with k=", Args[a+1], " did not finish"), sep="", collapse=""), stderr())
      not_finished = c(not_finished, Args[a+1])
      df = df[-nrow(df),]
   }
}

# checking whether any assembly finished #
if (all(is.na(df))) {
   stop("Not assemblies finished")
}

df[,"1/Nodes"] = 1/df[,"Nodes"]

# ranking N50, MaxContig, 1/Nodes
ranks = apply(df[,c(3,6,5)], 2, rank)
rownames(ranks) = df$Ksize

# summing ranks
ranks.s = apply(ranks, 1, sum)

# taking the one with highest ranks, if tie then chose by min no. nodes
m = max(ranks.s)
best = which(ranks.s == m)
if (length(best) != 1) {
   ks = names(best)
   best = ks[which.min(df[df$Ksize %in% ks,"Nodes"])]
   best_assembly = best
} else {
   best_assembly = names(best)
}

# plot
library(ggplot2)

col = rep("Not Best", nrow(df))
names(col) = df$Ksize
col[best_assembly] = "Best"
df$Status = col

df.melt = melt(df[,c(1,2,3,4,5)], id=c("Ksize"))
df.melt$Status = rep(col, 4)

p = ggplot(df.melt, aes(x=Ksize, y=value, colour=Status)) + facet_wrap(~variable, nrow=2, ncol=2, scales="free") + geom_point() 
pdf(file="velvet_parse.pdf", width=9, height=7)
print(p)
dev.off()

# write out
output = c(paste(c("Best assembly is: ", best_assembly), sep="", collapse=""),
           paste(c("Rank of assemblies (best-worst) (ties not properly sorted):", names(ranks.s[order(ranks.s, decreasing=TRUE)]), not_finished), sep=" ", collapse=" "))
write(output, file="velvet_parse.txt")

write.table(df, file="velvet_parse.tab", quote=FALSE, sep="\t", row.names=FALSE)


cat(paste(c("Best assembly is: ", best_assembly, "\n"), sep="", collapse=""))
cat(paste(c("Rank of assemblies (best-worst) (ties not properly sorted):", names(ranks.s[order(ranks.s, decreasing=TRUE)]), not_finished, "\n"), sep=" ", collapse=" "))



