source("R/utils.R")
source("R/mcmc.R")
source("R/functions.R")
load_pkgs()
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
theme_set(theme_bw())
library(targets)
library(phytools)
library(diversitree)


# FIGURE 1 ----------------------------------------------------------------
#simmap simulations
tar_load(ag_model_pcsc)
tar_load(treeblock)
tar_load(ag_compdata_tb)
phy<-treeblock[[1]]
data<-ag_compdata_tb$data
model<-ag_model_pcsc$solution

set.seed(1)
sims<-list()
sums<-list()
counts<-matrix(ncol=3, nrow = 0)
i<-1
while(i<101){
  sims[[i]]<-makeSimmap(tree = treeblock[[i]], data = data, model = model, rate.cat = 1, nSim = 10)
  sims[[i]]<-lapply(sims[[i]],mergeMappedStates,c(1:4),"ag0")
  sims[[i]]<-lapply(sims[[i]],mergeMappedStates,c(5:8),"ag1")
  class(sims[[i]])<-c("multiSimmap","multiPhylo")
  sums[[i]]<-summary(sims[[i]])
  counts<-rbind(counts,sums[[i]]$count)
  print(i)
  i=i+1
}

#gain and loss CIs
gain.ci<-quantile(counts[,2], c(0.025,0.975))
loss.ci<-quantile(counts[,3], c(0.025,0.975))

#finding AG nodes
nodeProbs<-as.data.frame(sums[[1]]$ace)
nodeProbs$ag<-as.factor(round(nodeProbs[,2], digits = 0))
allAgNodes<-rownames(subset(nodeProbs, ag==1))

#plotting posterior probability
obj1<-densityMap(sims[[1]], plot=FALSE)
n<-length(obj1$cols)
obj1$cols[1:n]<-colorRampPalette(c("grey60","firebrick"), space="Lab")(n)

plot(obj1,type="fan",ftype="off", lwd=2)
nodelabels(node = allAgNodes, pch = 21, col="firebrick4", bg="firebrick", cex=0.5, lwd=2)

#why cant i load this directly?
data<-read.csv(file.choose())
data2<-data.frame(care=data$care, spawning=data$spawning, row.names = data$species)

cols1<-list(care=c("#c7e9c0","#006d2c"), spawning=c("#bdd7e7","#2171b5"))
labs1<-c("Parental Care","Spawning Mode")
str1<-list(care=c("No Male Care","Male Care"), spawning=c("Pair Spawning", "Group Spawning"))

trait.plot(treeblock[[1]], data2, cols1, cex.lab = 0.2,lab = labs1, str = str1, cex.legend = 0.5)

# FIGURE 2 ----------------------------------------------------------------
tar_load(contr_long_ag_mcmc0)
tar_load(contr_long_ag_mcmc_tb)
tar_load(contr_long_ag_priorsamp)
ag_contr_gainloss <- (purrr::map_dfr(list(fishphylo=contr_long_ag_mcmc0,
                                          treeblock=contr_long_ag_mcmc_tb,
                                          prior = contr_long_ag_priorsamp),
                                     filter, rate != "netgain",
                                     .id = "phylo")
                      %>% mutate(across(phylo, factor,
                                        levels=rev(c("prior","treeblock","fishphylo"))))
                      ## but we only care about treeblock for the paper
                      %>% filter(phylo == "treeblock")
                      %>% filter(contrast != "intercept")
                      
)

ag_contr_gainloss$rate<-factor(ag_contr_gainloss$rate, levels = c("loss","gain"))

gg_sum_nice <- ggplot(ag_contr_gainloss, aes(x = exp(value), y = contrast, colour = rate)) +
  geom_violin(aes(fill = rate), alpha=0.6) +
  stat_summary(fun.data = "median_hilow",
               geom = "errorbar",
               aes(group=rate),
               ## width by trial and error; not sure what determines this?
               position = position_dodge(width=0.875),
               colour = "black") +
  stat_summary(fun = median,
               geom = "point", aes(group=rate),
               ## width by trial and error; not sure what determines this?
               position = position_dodge(width=0.875),
               colour = "black",
               pch = 3,
               size = 2) + 
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  zmargin +
  scale_colour_manual(name="", labels = c("Gain","Loss"), values = c("firebrick","gray70"), limits = c("gain", "loss"))+
  scale_fill_manual(name="", labels = c("Gain","Loss"), values = c("firebrick","gray70"), limits = c("gain", "loss"))+
  scale_y_discrete(breaks=c("pcxsc","sc","pc"), 
                   labels=c("Interaction", "Spawning Mode", "Parental Care"),
                   limits=c("pcxsc","sc","pc"))+
  labs(x="Proportional Difference in Rates", y="")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12))
print(gg_sum_nice)
