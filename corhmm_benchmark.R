load("cache/MK_3state_simple.rda")
## will be a benchmark of corHMM(..., p = ...)  vs running with arg list
nllfun(p)
fun2 <- function(p) {
  capture.output(cc <- corHMM(p = exp(p), phy = phy, data = data, rate.cat = 1, rate.mat = StateMatA_constrained))
  -(cc$loglik)
}
fun2(p)

library(rbenchmark)
b1 <- benchmark(nllfun(p), fun2(p), columns=c("test", "replications", "elapsed", "relative"))
print(b1)
## definitely worth returning arg.list: about 12x faster function evaluation
