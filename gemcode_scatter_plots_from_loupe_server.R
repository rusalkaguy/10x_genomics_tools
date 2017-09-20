#
# pull data (CSV) for structure plots from 10X Genomics Loupe server and graph
# for a list of samples
#

library(data.table)  # pull csv from URL

#
# loupe constants
#
loupeServerURL = "http://bmilinux.genome.uab.edu:3000/loupe/view/"
csvURL = "/loupe-sv-barcode-matrix.csv?"
axisLabels = c("chr1:161,496,873","chr1:161,845,411")
axisLabel = paste(axisLabels, collapse="-")
regionURL = "x=chr1%2B161496873-chr1%2B161845411&y=chr1%2B161496873-chr1%2B161845411"
csvOptionsURL = "&downsample=1&nr=yes&ns=yes&low_quality=no&restrict_haplotype_x=&restrict_haplotype_y="

# scraped from loupe HTML - see "Loupe Color Scales.xlsx" in box.com
colRamp29 = c("#FFFFFF"	,"#FFF0F0"	,"#FFE680"	,"#FFD240"	,"#FFC000"	,"#FFB000"	,"#FFA000"	,"#FF9000"	,"#FF8000"	,"#FF7000"	,"#FF6000"	,"#FF5000"	,"#FF4000"	,"#FF2000"	,"#FF0000"	,"#F00000"	,"#E00000"	,"#D20000"	,"#C40000"	,"#AE0000"	,"#A00000"	,"#900000"	,"#800000"	,"#700000"	,"#600000"	,"#500000"	,"#400000")

#
# index of our server's data
# created by hand, no way to query this from the server :-(
#
sampleKeys = c(
  "A_wt-hg38"            = "U0xFNDIwNV93dC1oZzM4LmxvdXBl"
  ,"B_2bv-hg38"          = "U0xFODA5NF8yYnYtaGczOC5sb3VwZQ=="
  ,"C_2bv-hg38"          = "U0xFNDA5MF8yYnYtaGczOC5sb3VwZQ=="
)

# create data cache
dataCacheFile = "loupe.datacache.RData"
if(file.exists(dataCacheFile)) { load(dataCacheFile, verbose=T) }
if(!exists("dataCache")){dataCache = list()}


# 
# set up images per page & pdf file
#
mfrow = c(3,2)
pdf(file=paste0("loupe.fcgr_cluster.heatmap.all.",paste(mfrow,collapse="x"),".pdf"), width=8.5,height=11)
opar=par(mfrow=mfrow)
for(sampleName in names(sampleKeys)) {
  # debug
  # sampleName="SLE1002_2bv-hg38"
  cat("#### ",sampleName, "####\n")
  url = paste0(loupeServerURL, sampleKeys[sampleName], csvURL, regionURL, csvOptionsURL)

  # lazy load
  mydat = dataCache[[url]]
  if(is.null(mydat)) {
    cat("Loading [",sampleName,"]: \n")
    # pull from internet, parse and cache
    dataCache[[url]] = fread(url, sep=',')
    mydat = dataCache[[url]]
  }
    
  (sums=summary(as.vector(as.matrix(mydat))))

  image(as.matrix(mydat)
        , col=colRamp29, main=sampleName
        #, xlab=axisLabel, ylab=axisLabel
        , xaxt="n", yaxt="n", pty="m"
  )
  segments(0,0,1,1, col="black")  # diagonal
  box() # bounding box
  
  #axis(1,at=c(0,1), labels=axisLabels)
  #axis(2,at=c(0,1), labels=axisLabels)
  axis(1,at=c(0.5), labels=axisLabel)
  axis(2,at=c(0.5), labels=axisLabel)
  
  #
  # heatmap doesn't work with mfrow()
  #
  # h=heatmap(as.matrix(mydat)
  #             , Rowv = NA, Colv=NA, labCol=NA, labRow=NA
  #             , main=sampleName, xlab=axisLabel, ylab=axisLabel, margins=c(0,0)
  #             , scale="none"  # prevent row or col centering!!!
  #             , col=colRamp29) 
  # 
  

}
par(opar)
dev.off()

#
# save data for next time
#
save(file=dataCacheFile, dataCache)
