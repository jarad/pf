source("pf_plyr_functions.r")
require(plyr)

mydata = expand.grid(n = c(100,1000,10000,20000), filt=c("KD","RM"), resamp="systematic", prior="lognormal", mod.data="orig", mod.fit="orig", progress=TRUE, stringsAsFactors=FALSE)
mydata = data.frame(mydata, seed = 61:64)

require(doMC)
registerDoMC()

m_ply(.data = mydata, .fun = pf, .parallel = TRUE)
m_ply(.data = mydata[,-((dim(mydata)[2]-1):(dim(mydata)[2]))], .fun = pf.quant, .parallel = TRUE)
