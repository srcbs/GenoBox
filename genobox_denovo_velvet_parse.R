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
   
   # N50
   N50 = calcN50.f(contig_lengths)
   
   return (c(contig_sum, N50))
}


#Args = c("", "", "test_33/stats.txt", 33, "test_37/stats.txt", 37, "test_41/stats.txt", 41, "test_45/stats.txt", 45, "test_49/stats.txt", 49)
#Args = c("/tools/bin/R-2.12", "--vanilla", "Kleb_33/stats.txt", "33", "Kleb_37/stats.txt", "37", "Kleb_41/stats.txt", "41", "Kleb_45/stats.txt", "45", "Kleb_49/stats.txt", "49", "Kleb_53/stats.txt", "53", "Kleb_57/stats.txt", "57", "Kleb_61/stats.txt", "61", "Kleb_65/stats.txt", "65", "Kleb_69/stats.txt", "69", "Kleb_73/stats.txt", "73")

# parse inputs to matrix
row = 0
df = as.data.frame(matrix(NA, ncol=4, nrow=(length(Args)-2)/2))
colnames(df) = c("Ksize", "Sum", "N50", "Sum_x_N50")
for (a in seq(3,length(Args),2)) {
   row = row + 1
   df[row,1] = Args[a+1]
   df[row,c(2,3)] = parse.f(Args[a], as.integer(Args[a+1]))
}

# multiply the Sum and N50 to get best assembly:
df[,4] = df[,2]*df[,3]

rank_assembly = df[order(df[,4], decreasing=TRUE),1]

# plot
library(ggplot2)
p = ggplot(df, aes(x=Sum, y=N50,label=Ksize)) + geom_text()
pdf(file="velvet_parse.pdf", width=7, height=7)
print(p)
dev.off()

# write out
output = c(paste(c("Best assembly is: ", rank_assembly[1]), sep="", collapse=""),
           paste(c("Rank of assemblies (best-worst):", rank_assembly), sep=" ", collapse=" "))
write(output, file="velvet_parse.txt")

write.table(df, file="velvet_parse.tab", quote=FALSE, sep="\t", row.names=FALSE)


cat(paste(c("Best assembly is: ", rank_assembly[1], "\n"), sep="", collapse=""))
cat(paste(c("Rank of assemblies (best-worst):", rank_assembly, "\n"), sep=" ", collapse=" "))



