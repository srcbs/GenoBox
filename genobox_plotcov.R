# plot coverage from bedtools as histogram and as cumulative plot. Excludes 0 covered positions

# input are <infile> <title> <outfile>
cat("-- reading arguments\n", sep = "")
Args <- commandArgs();

cov = read.csv(Args[3], sep="\t", as.is=TRUE, header=FALSE)
cov = cov[cov[,1] == "genome",]
cov.df = cov[-1,]

library(ggplot2)

# calc upper depth limit for x-axis
upper = max((which((cov.df[,3]/cov.df[,4])*100 > 0.1)))

# create standard histogram
p1 = ggplot(cov.df, aes(V2, V3)) + geom_bar(stat="identity") + xlim(1,upper) + opts(title=Args[4]) + labs(x="Depth", y="Observations")

# create cumulative plot

# fraction covered with more than 10
csums = sum(cov.df[10:nrow(cov.df),5])
t = rev(cumsum(rev(cov.df$V5)))*100
df = data.frame(Depth=1:length(t), Percentage=t)
p2 = ggplot(df, aes(Depth, Percentage)) + geom_line() + xlab(">= X depth") + ylab("% genome covered") + scale_y_continuous(limits=c(0,100)) + xlim(1,2*upper) + opts(title=Args[4])

# plot
pdf(file=Args[5], width=10, height=7)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))
vplayout <- function(x,y)
viewport(layout.pos.row=x, layout.pos.col=y)
print(p1, vp=vplayout(1,1))
print(p2, vp=vplayout(1,2))
dev.off()

# print out fraction covered at 10X or more
cat(paste(c("Covered at 10X or more: ", signif(csums*100, digits=3), "\n"), sep="", collapse=""))
