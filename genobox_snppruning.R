# read commandline arguments
cat("-- reading arguments\n", sep = "")
Args <- commandArgs();

# commandline arguments are <N bases to discard SNPs within> <inputfile> <outputfile>
#Args = c("","",5, "Aborigine_trimmed_hg19.var.snps.altfilt.af.vcf", "Aborigine_trimmed_hg19.var.snps.altfilt.af.pruned_5nt.vcf")

vcf = read.csv(Args[4], sep="\t", comment.char="#", header=FALSE, as.is=TRUE)

# calc difference to next snp
d_list = list()
vcf_list = list()
for (c in unique(vcf[,1])) {
   vcf_list[[c]] = vcf[vcf[,1] == c,]
   d_list[[c]] = diff(vcf_list[[c]][,2])
}
names(vcf_list) = unique(vcf[,1])
names(d_list) = unique(vcf[,1])

# remove from sets
vcf_filtered_list = lapply(names(d_list), function(i) {
   d = d_list[[i]]
   set = vcf_list[[i]]
   rms = which(d <= as.numeric(Args[3]))
   index = c()
   if (length(rms) == 0) {
      set
   } else {
      for (r in rms) {
         index = c(index, r:(r+1))
      }
      index = unique(index)
      set_filtered = set[-index,]
      set_filtered
   }
})
names(vcf_filtered_list) = names(d_list)

# create output-table
for (i in 1:length(vcf_filtered_list)) {
   if (i == 1) {
      out.df = vcf_filtered_list[[i]][,1:10]
   } else {
      out.df = rbind(out.df, vcf_filtered_list[[i]][,1:10])
   }
}

# write out
filtered = nrow(vcf) - nrow(out.df)
out.txt = paste(c("Total ", nrow(vcf), ", Written ", nrow(out.df), ", Filtered ", filtered, "\n"), sep="", collapse="")
cat(out.txt)
write.table(out.df, file=Args[5], row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
