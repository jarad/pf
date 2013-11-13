source("pf_plyr_functions.r")
require(plyr)

pf(100, "RM", "systematic", "uniform", "orig", "orig")
pf.quant(100, "RM", "systematic", "uniform", "orig", "orig")

mydata = expand.grid(n = c(100,1000,10000,20000), filt="KD", resamp="systematic", prior="uniform", mod.data="orig", mod.fit="orig", progress=FALSE, stringsAsFactors=FALSE)

# Apply pf to combination of pfs
#data1 = expand.grid(mod.dat = "orig", mod.fit = "orig", n = c(100, 1000, 10000, 20000, 40000), filt = c("BF","APF","KD"), resamp = "systematic", prior = "uniform", progress=FALSE, stringsAsFactors=FALSE)
#data2 = data.frame(mod.dat = "orig", mod.fit = "orig", n = c(10000, 20000), filt = "KD", resamp = "systematic", prior = "semi-uniform", progress=FALSE, stringsAsFactors=FALSE)
#data3 = expand.grid(mod.dat = "orig", mod.fit = "orig", n = c(100, 1000, 10000, 20000, 40000), filt = "KD", resamp = c("multinomial","residual","stratified","systematic"), prior = "lognormal", progress=FALSE, stringsAsFactors=FALSE)
#data4 = expand.grid(mod.dat = "ext", mod.fit = c("ext","orig"), n = c(100, 1000, 10000, 20000, 40000), filt = "KD", resamp = "stratified", prior = "lognormal", prior.samp = "lognormal", progress=FALSE, stringsAsFactors=FALSE)
#seed = c(rep(61:65,3), 63, 64, rep(61:65, 4))
#mydata = data.frame(rbind(data1,data2,data3),seed)

require(doMC)
registerDoMC()

m_ply(.data = mydata, .fun = pf, .parallel = TRUE)
m_ply(.data = mydata[,-dim(mydata)[2]], .fun = pf.quant, .parallel = TRUE)