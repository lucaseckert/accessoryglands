load("cache/MK_3state_simple.rda")
library(corHMM)
## will be a benchmark of corHMM(..., p = ...)  vs running with arg list
nllfun(p)  ## BMB version
fun2 <- function(log_p) {
  capture.output(cc <- corHMM(p = exp(log_p), phy = phy, data = data, rate.cat = 1, rate.mat = StateMatA_constrained))
  -(cc$loglik)
}
fun2(p)

library(rbenchmark)
b1 <- benchmark(nllfun(p), fun2(p), columns=c("test", "replications", "elapsed", "relative"))
print(b1)
## definitely worth returning arg.list: about 12x faster function evaluation
