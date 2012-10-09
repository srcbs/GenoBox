# read commandline arguments
cat("-- reading arguments\n", sep = "")
Args <- commandArgs();

# commandline arguments are <Map_file> <Rmdup_file> <Name>
#Args = c("","", "Athabascan_1_1/Athabascan_1_1_Map.txt", "Athabascan_1_1/Athabascan_1_1_Rm.txt", "Athabascan_1_1_sub")
#Args = c("","", "UM29.13082012/UM29.13082012_Map.txt", "UM29.13082012/UM29.13082012_Rm.txt", "UM29.13082012")

Map = read.table(Args[3], header=FALSE)
Rm = read.table(Args[4], header=FALSE)
lim=max(c(max(Map[,3]), max(Rm[,3])))
pdf(file=paste(c(Args[5], ".saturation.pdf"), sep="", collapse="") , width=7, height=7)
plot(Map[,3], Rm[,3], main=Args[5], xlab="Mapped reads", ylab="Rmduped reads", pch=16, xlim=c(0,lim), ylim=c(0,lim))
abline(a=0,b=1, col="grey")
dev.off()

v = Rm[,3]/Map[,3]
output = c(Args[5], mean(v), sd(v))
write.table(t(as.matrix(output)), file=paste(c(Args[5], ".stats.txt"), sep="", collapse=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

mean(v)
sd(v)