
# implement Moab class (currently done for alignment) # ALIGNMENT, BAMPROCESS, BAMSTATS, GENOTYPING, VCFFILTER, DBSNP

# implement libfile as class
# move library stuff from genobox_modules so that it is inside the class

genobox.py bcf2ref --bcf abcalls.all.bcf --genome /panvol1/simon/projects/aborigine_stinus/hg19_flat/build37_rCRS.genome --ex /panvol1/simon/databases/hs_ref37_rCRS/gi2number.build37_rCRS --dbsnp /panvol1/simon/databases/dbsnp/dbsnp132_hg19.vcf.gz --rmsk /panvol1/simon/databases/hs_ref37_rCRS/rmsk/rmsk_build37_rCRS.number.sort.genome --indels indels_for_filtering.number.vcf