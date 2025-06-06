source("R/utils.R")
source("R/mcmc.R")
source("R/functions.R")
source("R/pkg_install.R")
load_pkgs()
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
theme_set(theme_bw())
library(targets)
library(phytools)
library(diversitree)

## FIXME: reconcile themes/DRY?

## save sim 2, most "typical"
sims_save <- c(2)
### code for simulations 
### I just saved the output as rds since I'm not sure how to work the targets
# tar_load(ag_model_pcsc)
# tar_load(treeblock)
# tar_load(ag_compdata_tb)
# data<-ag_compdata_tb$data
# model<-ag_model_pcsc$solution

# set.seed(1)
# sims<-list()
# sums<-list()
# counts<-matrix(ncol=3, nrow = 0)
# i<-1
# while(i<101){
#   sims[[i]]<-makeSimmap(tree = treeblock[[i]], data = data, model = model, rate.cat = 1, nSim = 100)
#   sims[[i]]<-lapply(sims[[i]],mergeMappedStates,c(1:4),"ag0")
#   sims[[i]]<-lapply(sims[[i]],mergeMappedStates,c(5:8),"ag1")
#   class(sims[[i]])<-c("multiSimmap","multiPhylo")
#   sums[[i]]<-summary(sims[[i]])
#   counts<-rbind(counts,sums[[i]]$count)
#   print(i)
#   i=i+1
# }

#simmap simulations
tar_load(ag_model_pcsc)
tar_load(treeblock)
tar_load(ag_compdata_tb)
phy <- treeblock[[1]]
data <- ag_compdata_tb$data
model <- ag_model_pcsc$solution

set.seed(1)
sims <- list()
sums <- list()
counts <- list()
i <- 1

sim0 <- makeSimmap(tree = phy,
                   data = data, model = model,
                   rate.cat = 1, nSim = 1)[[1]]
statenames <- rownames(sim0$Q)
nsims <- 100

if (file.exists("simmap.rda")) {
    load("simmap.rda")
} else {
    pb <- txtProgressBar(max = nsims, style = 3)
    for (i in 1:length(treeblock)) {
        setTxtProgressBar(pb, i)
        sims[[i]] <- makeSimmap(tree = treeblock[[i]],
                            data = data, model = model,
                            rate.cat = 1, nSim = nsims)
        sims[[i]] <- lapply(sims[[i]], mergeMappedStates, statenames[1:4], "ag0")
        sims[[i]] <- lapply(sims[[i]], mergeMappedStates, statenames[5:8], "ag1")
        ## restore state (dropped by lapply())
        class(sims[[i]]) <- c("multiSimmap","multiPhylo")
        sums[[i]] <- summary(sims[[i]])
        dim(sums[[i]]$count)
        counts[[i]] <- sums[[i]]$count
    }
    close(pb)
    counts <- do.call("rbind", counts)
    ## simplify:
    sims <- sims[sims_save]
    sums <- sums[sims_save]
    save(counts, sums, sims, file = "simmap.rda")
}

#gain and loss CIs
gain.ci <- quantile(counts[,"ag0,ag1"], c(0.025,0.975))
loss.ci <- quantile(counts[,"ag1,ag0"], c(0.025,0.975))

## BMB: where are these used
?
#finding AG nodes

which_sim <- 1 ## this is indexed based on *which sims were saved* (not original 1:100)
nodeProbs <- as.data.frame(sums[[which_sim]]$ace) 
nodeProbs$ag <- as.factor(round(nodeProbs[,2], digits = 0))
allAgNodes <- rownames(subset(nodeProbs, ag==1))

#plotting posterior probability
obj1 <- densityMap(sims[[which_sim]], plot=FALSE)
n <- length(obj1$cols)
obj1$cols[1:n] <- colorRampPalette(c("grey60","firebrick"), space="Lab")(n)

agNodes<-c(778,784,794,798,808,817,876,923,1018,1020,1103,"Hoplosternum_littorale",
           "Ompok_siluroides","Pangasius_pangasius","Auchenipterus_nuchalis",
           "Lepidogalaxias_salamandroides","Cheilinus_undulatus","Radulinopsis_taranetzi")

plot(obj1,type="fan",ftype="off", lwd=2)
nodelabels(node = agNodes, pch = 21, col="firebrick4", bg="firebrick", cex=1.125, lwd=2)
tiplabels(tip = agNodes, pch = 21, col="black", bg="firebrick", cex=1, lwd=2)



## why cant i load this directly from targets?
data <- read.csv("data/binaryTraitData.csv")
data2 <- data[c("care", "spawning")]
rownames(data2) <- data$species

cols1 <- list(care=c("#c7e9c0","#006d2c"), spawning=c("#bdd7e7","#2171b5"))
labs1 <- c("Male Care","Spawning Mode")
str1 <- list(care=c("No Male Care","Male Care"), spawning=c("Pair Spawning", "Group Spawning"))

do_png <- TRUE
tar_load(treeblock)
totpix <- 4800
res <- 600
tfun <- function() {
    trait.plot(treeblock[[2]], data2, cols1, cex.lab = 0.2,
               lab = labs1, str = str1, cex.legend = 0.5)
}

png("ridiculous_tree.png", width = totpix, height = totpix, res = res,
    type = "cairo-png")
tfun()
dev.off()
file.info("ridiculous_tree.png")

svg("ridiculous_tree.svg", width = 8, height = 8)
tfun()
dev.off()
file.info("ridiculous_tree.svg")

#distribution of transitions
transitions<-data.frame(counts[,-1])
transitions$sim<-1:10000
transitions<-transitions %>% rename(
  gains=ag0.ag1,
  losses=ag1.ag0
)
trans<-transitions %>% gather(rate, n, gains,losses)

ggplot(trans, aes(x=n,fill=rate,))+
  geom_histogram(position = "identity", alpha=0.65, color="black", binwidth = 1)+
  scale_fill_manual(name="", labels = c("Gains","Losses"), values = c("firebrick","gray70"), limits = c("gains", "losses"))+
  labs(x="Transitions", y="Frequency")+
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size=12))

# FIGURE 2 ----------------------------------------------------------------
tar_load(contr_long_ag_mcmc0)
tar_load(contr_long_ag_mcmc_tb)
tar_load(contr_long_ag_priorsamp)
ag_contr_gainloss  <-  (purrr::map_dfr(list(fishphylo=contr_long_ag_mcmc0,
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

ag_contr_gainloss$rate <- factor(ag_contr_gainloss$rate, levels = c("loss","gain"))

cvec <- c("firebrick","gray70")
nvec <- c("Gain", "Loss")
ylab_pos <- c(2.2, 1.7) ## locations for gain/loss labels
vw <- 0.5 ## violin width
pdw <- 0.49  ## dodging (trial and error; depends on violin width)
gg_sum_nice <- ggplot(ag_contr_gainloss, aes(x = exp(value), y = contrast, colour = rate)) +
  geom_violin(aes(fill = rate), alpha=0.6, width = vw) +
  stat_summary(fun.data = "median_hilow",
               geom = "errorbar",
               width  = 0.1,
               aes(group=rate),
               ## width by trial and error; not sure what determines this?
               position = position_dodge(width=pdw),
               colour = "black") +
    stat_summary(fun = median,
                 geom = "point", aes(group=rate),
                 ## width by trial and error; not sure what determines this?
                 position = position_dodge(width=pdw),
                 colour = "black",
                 pch = 3,
                 size = 2) + 
    geom_vline(xintercept = 1, lty = 2) +
    scale_x_log10(labels = function(x) {
        ## ugly: don't want trailing zeros
        trimws(gsub("\\.00$", "", format(x, scientific = FALSE)))
    }) + 
    zmargin +
  scale_colour_manual(name="", labels = c("Gain","Loss"), values = cvec, limits = c("gain", "loss"))+
  scale_fill_manual(name="", labels = c("Gain","Loss"), values = cvec, limits = c("gain", "loss"))+
  annotate("text", label = nvec, col = cvec, x = 15, y = ylab_pos, size = 8) +
  scale_y_discrete(breaks=c("pcxsc","sc","pc"), 
                   labels=c("Interaction", "Group Spawning", "Male Care"),
                   limits=c("pcxsc","sc","pc"))+
  labs(x="Proportional Difference in Rates", y="")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 16),
        legend.position = "none")
print(gg_sum_nice)

ggsave("fig2.png", width = 7, height  = 5)

## 
system("eog fig2.png & ")


### Figures for Presentations ###
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
                      %>% filter(rate == "loss")
                      
)
gg_sc <- ggplot(ag_contr_gainloss, aes(x = exp(value), y = contrast, colour = rate)) +
  geom_violin(aes(fill = rate), alpha=0.6) +
  stat_summary(fun.data = "median_hilow",
               geom = "errorbar",
               ## width by trial and error; not sure what determines this?
               position = position_dodge(width=0.875),
               colour = "black") +
  stat_summary(fun = median,
               geom = "crossbar", aes(group=rate),
               ## width by trial and error; not sure what determines this?
               position = position_dodge(width=0.875),
               colour = "black",
               size = 1) +
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
        axis.text = element_text(size = 24, color = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.y = element_text(size = 24, color = "black"))

print(gg_sc)

##Proportion Figures##
data<-read.csv("data/binaryTraitData.csv")

care<-subset(data, care==1)
noCare<-subset(data, care==0)

AgWithCare<-sum(care$ag)/length(care$ag)
noAgWithCare<-sum(noCare$ag)/length(noCare$ag)

groups<-subset(data, spawning==1)
pairs<-subset(data, spawning==0)

AgGroups<-sum(groups$ag)/length(groups$ag)
AgPairs<-sum(pairs$ag)/length(pairs$ag)

careProps<-c(AgWithCare, noAgWithCare)
careStates<-c("Male Care", "No Male Care")
careData<-data.frame(careStates, careProps)

spawnProps<-c(AgGroups, AgPairs)
spawnStates<-c("Groups", "Pairs")
spawnData<-data.frame(spawnStates, spawnProps)

ggplot(careData, aes(x=careStates, y=careProps, fill=careStates))+
  geom_bar(stat = "identity", color="black")+
  scale_fill_manual(name="", labels = c("care","none"), 
                    values = c("gray30","gray70"), 
                    limits = c("Male Care", "No Male Care"))+
  labs(x="",y="Proportion with Accessory Glands")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.position = "none")

ggplot(spawnData, aes(x=spawnStates, y=spawnProps, fill=spawnStates))+
  geom_bar(stat = "identity", color="black")+
  scale_fill_manual(name="", labels = c("groups","pairs"), 
                    values = c("gray30","gray70"), 
                    limits = c("Groups", "Pairs"))+
  labs(x="",y="Proportion with Accessory Glands")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.position = "none")
