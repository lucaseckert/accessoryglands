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
set.seed(1)
sm1 <- with(ag_model_pcsc,
            makeSimmap(phy, data, solution, rate.cat, nSim = 100))

#finding AG nodes
smSum<-summary(sm1)
nodeProbs<-as.data.frame(smSum$ace)
nodeProbs$ag1<-rowSums(nodeProbs[,5:8])
nodeProbs$ag<-as.factor(round(nodeProbs$ag1, digits = 0))
allAgNodes<-rownames(subset(nodeProbs, ag==1))

#collapsing to AG
agSims<-lapply(sm1,mergeMappedStates,c(1:4),"ag0")
agSims<-lapply(agSims,mergeMappedStates,c(5:8),"ag1")
class(agSims)<-c("multiSimmap","multiPhylo")

#plotting posterior probability
obj1<-densityMap(agSims, plot=FALSE)
n<-length(obj1$cols)
obj1$cols[1:n]<-colorRampPalette(c("grey60","firebrick"), space="Lab")(n)

#The node/tips where AGs first evolved: Lepidogalaxias salamandroides, 767, Cheilinus undulatus, 737, 615, 533, Hoplosternum littorale, 850, 845, Rita rita, 839, 836, Pangasius pangasius, 830
agNodes<-c(767,737,615,533,850,845,839,836,830)
agTips<-c("Lepidogalaxias_salamandroides","Cheilinus_undulatus","Hoplosternum_littorale","Rita_rita","Pangasius_pangasius")

plot(obj1,type="fan",ftype="off", lwd=2)
nodelabels(node = agNodes, pch = 21, col="firebrick4", bg="firebrick", cex=1.5, lwd=2)
tiplabels(tip = agTips, pch = 21, col="firebrick4", bg="firebrick", cex=1.5, lwd=2)
#tip labels not working for some reason

#why cant i load this directly?
data<-read.csv(file.choose())
data2<-data.frame(care=data$care, spawning=data$spawning, row.names = data$species)

cols1<-list(care=c("#c7e9c0","#006d2c"), spawning=c("#bdd7e7","#2171b5"))
labs1<-c("Parental Care","Spawning Mode")
str1<-list(care=c("No Male Care","Male Care"), spawning=c("Pair Spawning", "Group Spawning"))

trait.plot(ag_model_pcsc$phy, data2, cols1, cex.lab = 0.2,lab = labs1, str = str1, cex.legend = 0.5)


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

# SIMMAPs -----------------------------------------------------------------
tar_load(ag_model_pcsc)
tar_load(treeblock)
set.seed(1)
root<-c(1,0,0,0,0,0,0,0)

sm1 <- makeSimmap(treeblock[[1]], ag_model_pcsc$data, ag_model_pcsc$solution, ag_model_pcsc$rate.cat, root.p = "yang", nSim = 1)

