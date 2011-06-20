library(cqnpaper)
samp <- c("NA11918", "NA12761")
data(montgomery.uCounts)
montgomery.subset <- montgomery.uCounts[, samp]
save(montgomery.subset, file = "montgomery.subset.rda")
data(sizeFactors.montgomery)
sizeFactors.subset <- sizeFactors.montgomery[samp]
save(sizeFactors.subset, file = "sizeFactors.subset.rda")
data(uCovar)
uCovar <- uCovar[,1:2]
save(uCovar, file = "uCovar.rda")

samp <- c("NA06985", "NA06994", "NA07037", "NA10847", "NA11920",
          "NA11918", "NA11931", "NA12003", "NA12006", "NA12287")
montgomery.subset <- montgomery.uCounts[, samp]
wh <- which(rowSums(montgomery.subset) > 0)
montgomery.subset <- montgomery.subset[wh,]
save(montgomery.subset, file = "montgomery.subset.rda")
data(sizeFactors.montgomery)
sizeFactors.subset <- sizeFactors.montgomery[samp]
save(sizeFactors.subset, file = "sizeFactors.subset.rda")
data(uCovar)
uCovar <- uCovar[wh,1:2]
save(uCovar, file = "uCovar.rda")
