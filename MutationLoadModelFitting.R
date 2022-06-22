## Model fitting for Mutation Load projec

#  Loading Packages 
require("matrixcalc")
require("coxme") # This package contains the lmekin function for fitting a mixed model
require("gdata") # Contains a function to fill in the upper/lower traingle of a symmetrical matrix
require("lqmm") # This package contains a function to make a matrix positive-semidefinite
require("lme4")
require("car")
require("ggplot2")
require("effects")
require("emmeans")
require("sjPlot")
require("sjmisc")


# Loading functions
# This is a function that will add a small value to the diagonal to make the kinship matrix positive semi-definite
normalize_kinmat <- function(kinmat){
  #normalize kinship so that Kij \in [0,1]
  tmp=kinmat - min(kinmat)
  tmp=tmp/max(tmp)
  tmp[1:9,1:9]
  #fix eigenvalues to positive
  diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
  tmp[1:9,1:9]
  return(tmp)
}

##  This is a function to perform an likelihood ratio test on the lmekin LMM output to extract a p-value
lmekin.anova <- function(m0,m1) {
  if((m1$method=="ML") & (m0$method=="ML")){
    k0=length(m0$coefficients$random)+length(m0$coefficients$fixed)
    k1=length(m1$coefficients$random)+length(m1$coefficients$fixed)
    k=k1-k0
    if(k>1){
      print("Are you sure these models are nested?")
    }
    D=2*(m1$loglik - m0$loglik)
    p=signif(1-pchisq(D,k), 2)
  }else{
    warning("Models need to be fit using ML")
    p=NA
  }
  return(p)
}

######################################################################################################################################################

############################################
######### Analysis with DSPR Data  #########
############################################

# Loading in phenotypic data with mutation load info
dspr.data = read.delim("/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Datasets/DSPR.Big.Data.Table.Feb.18.txt", header = TRUE)
dspr.data = dspr.data[order(dspr.data$DSPR.stock.number),]
dspr.data$DSPR.stock.number = paste("X", dspr.data$DSPR.stock.number, sep = "")

# Loading in fitness data with individual replicate measurements
dspr.fitness.data = read.csv("/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Datasets/transformed.full.data.oct.2019.csv")
dspr.fitness.data$DSPR.stock.number = paste("X", dspr.fitness.data$DSPR.stock.number, sep = "")

# Loading in genetic relatedness matrix matrix # This is a kinship matrix based on relatedness of haplotype blocks, not SNPs
dspr.GRM = read.delim("/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/KinshipMatricies/kinship.matrix.all.rils.txt", header = T, sep = "\t")
#names(dspr.GRM) = substring(names(dspr.GRM), 2,6)
rownames(dspr.GRM) = colnames(dspr.GRM)
#Preparing GRM
lowerTriangle(dspr.GRM) <- upperTriangle(dspr.GRM, byrow=TRUE) #This is a function that will convert the upper triangle to the lower triangle

# Normalize the GRM
# lhm.GRM.smoothed = as.matrix(cor.smooth(lhm.GRM)) #This is a function that will make the GRM positive definite
dspr.GRM = make.positive.definite(dspr.GRM) #This is a function that will make the GRM positive definite
# There is an additional RIL ID that doesn't exist in the pheno data that I remove here
dspr.GRM = dspr.GRM[rownames(dspr.GRM) %in% dspr.data$DSPR.stock.number, colnames(dspr.GRM) %in% dspr.data$DSPR.stock.number]

# Make a dummy matrix with each element equal to 1
dspr.GRM.null = matrix(1, ncol = length(rownames(dspr.GRM)), nrow = length(rownames(dspr.GRM)))
colnames(dspr.GRM.null) = rownames(dspr.GRM)
rownames(dspr.GRM.null) = rownames(dspr.GRM)

## Formatting phenotype data for model - formatting sex column
dspr.data.male.tmp = dspr.data[, -c(4:13)]
dspr.data.female.tmp = dspr.data[, -c(2,3,6:13)]
names(dspr.data.male.tmp)[2] = "fitness.cage"
names(dspr.data.male.tmp)[3] = "fitness.vial"
names(dspr.data.female.tmp)[2] = "fitness.cage"
names(dspr.data.female.tmp)[3] = "fitness.vial"
dspr.data.male.tmp$sex = "m"
dspr.data.female.tmp$sex = "f"
dspr.data.reformatted = rbind(dspr.data.male.tmp, dspr.data.female.tmp)

## Formatting phenotype data for model - formatting environment column
dspr.data.vial.tmp = dspr.data.reformatted[, -c(2)]
dspr.data.cage.tmp = dspr.data.reformatted[, -c(3)]
names(dspr.data.vial.tmp)[2] = "fitness"
names(dspr.data.cage.tmp)[2] = "fitness"
dspr.data.vial.tmp$mating.regime = "vial"
dspr.data.cage.tmp$mating.regime = "cage"
dspr.data.reformatted = rbind(dspr.data.vial.tmp, dspr.data.cage.tmp)

## Scale each trait by performing a Z transformation to each
dspr.data.reformatted$most.severe.variants.whole.genome.z.transformed = (dspr.data.reformatted$most.severe.variants.whole.genome - mean(dspr.data.reformatted$most.severe.variants.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$most.severe.variants.whole.genome,na.rm = TRUE)
dspr.data.reformatted$synonymous.variants.whole.genome.z.transformed = (dspr.data.reformatted$synonymous.variants.whole.genome - mean(dspr.data.reformatted$synonymous.variants.whole.genome, na.rm = TRUE)) / sd(dspr.data.reformatted$synonymous.variants.whole.genome, na.rm=TRUE)
dspr.data.reformatted$missense.variants.whole.genome.z.transformed = (dspr.data.reformatted$missense.variants.whole.genome - mean(dspr.data.reformatted$missense.variants.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$missense.variants.whole.genome,na.rm = TRUE)
dspr.data.reformatted$stop.gained.variants.whole.genome.z.transformed = (dspr.data.reformatted$stop.gained.variants.whole.genome - mean(dspr.data.reformatted$stop.gained.variants.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$stop.gained.variants.whole.genome,na.rm = TRUE)
dspr.data.reformatted$stop.lost.variants.whole.genome.z.transformed = (dspr.data.reformatted$stop.lost.variants.whole.genome - mean(dspr.data.reformatted$stop.lost.variants.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$stop.lost.variants.whole.genome,na.rm = TRUE)
dspr.data.reformatted$splice.disrupting.variants.whole.genome.z.transformed = (dspr.data.reformatted$splice.disrupting.variants.whole.genome - mean(dspr.data.reformatted$splice.disrupting.variants.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$splice.disrupting.variants.whole.genome,na.rm = TRUE)
dspr.data.reformatted$zhang.radical.score.polarity.and.volume.whole.genome.z.transformed = (dspr.data.reformatted$zhang.radical.score.polarity.and.volume.whole.genome - mean(dspr.data.reformatted$zhang.radical.score.polarity.and.volume.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$zhang.radical.score.polarity.and.volume.whole.genome,na.rm = TRUE)
dspr.data.reformatted$INDELs.whole.genome.z.transformed = (dspr.data.reformatted$INDELs.whole.genome - mean(dspr.data.reformatted$INDELs.whole.genome,na.rm = TRUE)) / sd(dspr.data.reformatted$INDELs.whole.genome,na.rm = TRUE)

dspr.data.reformatted$most.severe.variants.autosome.z.transformed = (dspr.data.reformatted$most.severe.variants.autosome - mean(dspr.data.reformatted$most.severe.variants.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$most.severe.variants.autosome,na.rm = TRUE)
dspr.data.reformatted$synonymous.variants.autosome.z.transformed = (dspr.data.reformatted$synonymous.variants.autosome - mean(dspr.data.reformatted$synonymous.variants.autosome, na.rm = TRUE)) / sd(dspr.data.reformatted$synonymous.variants.autosome, na.rm=TRUE)
dspr.data.reformatted$missense.variants.autosome.z.transformed = (dspr.data.reformatted$missense.variants.autosome - mean(dspr.data.reformatted$missense.variants.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$missense.variants.autosome,na.rm = TRUE)
dspr.data.reformatted$stop.gained.variants.autosome.z.transformed = (dspr.data.reformatted$stop.gained.variants.autosome - mean(dspr.data.reformatted$stop.gained.variants.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$stop.gained.variants.autosome,na.rm = TRUE)
dspr.data.reformatted$stop.lost.variants.autosome.z.transformed = (dspr.data.reformatted$stop.lost.variants.autosome - mean(dspr.data.reformatted$stop.lost.variants.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$stop.lost.variants.autosome,na.rm = TRUE)
dspr.data.reformatted$splice.disrupting.variants.autosome.z.transformed = (dspr.data.reformatted$splice.disrupting.variants.autosome - mean(dspr.data.reformatted$splice.disrupting.variants.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$splice.disrupting.variants.autosome,na.rm = TRUE)
dspr.data.reformatted$zhang.radical.score.polarity.and.volume.autosome.z.transformed = (dspr.data.reformatted$zhang.radical.score.polarity.and.volume.autosome - mean(dspr.data.reformatted$zhang.radical.score.polarity.and.volume.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$zhang.radical.score.polarity.and.volume.autosome,na.rm = TRUE)
dspr.data.reformatted$INDELs.autosome.z.transformed = (dspr.data.reformatted$INDELs.autosome - mean(dspr.data.reformatted$INDELs.autosome,na.rm = TRUE)) / sd(dspr.data.reformatted$INDELs.autosome,na.rm = TRUE)

dspr.data.reformatted$most.severe.variants.X.chromosome.z.transformed = (dspr.data.reformatted$most.severe.variants.X.chromosome - mean(dspr.data.reformatted$most.severe.variants.X.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$most.severe.variants.X.chromosome,na.rm = TRUE)
dspr.data.reformatted$synonymous.variants.X.chromosome.z.transformed = (dspr.data.reformatted$synonymous.variants.X.chromosome - mean(dspr.data.reformatted$synonymous.variants.X.chromosome, na.rm = TRUE)) / sd(dspr.data.reformatted$synonymous.variants.X.chromosome, na.rm=TRUE)
dspr.data.reformatted$missense.variants.X.chromosome.z.transformed = (dspr.data.reformatted$missense.variants.X.chromosome - mean(dspr.data.reformatted$missense.variants.X.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$missense.variants.X.chromosome,na.rm = TRUE)
dspr.data.reformatted$stop.gained.variants.X.chromosome.z.transformed = (dspr.data.reformatted$stop.gained.variants.X.chromosome - mean(dspr.data.reformatted$stop.gained.variants.X.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$stop.gained.variants.X.chromosome,na.rm = TRUE)
dspr.data.reformatted$stop.lost.variants.X.chromosome.z.transformed = (dspr.data.reformatted$stop.lost.variants.X.chromosome - mean(dspr.data.reformatted$stop.lost.variants.X.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$stop.lost.variants.X.chromosome,na.rm = TRUE)
dspr.data.reformatted$splice.disrupting.variants.X.chromosome.z.transformed = (dspr.data.reformatted$splice.disrupting.variants.X.chromosome - mean(dspr.data.reformatted$splice.disrupting.variants.X.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$splice.disrupting.variants.X.chromosome,na.rm = TRUE)
dspr.data.reformatted$zhang.radical.score.polarity.and.volume.X.chromosome.z.transformed = (dspr.data.reformatted$zhang.radical.score.polarity.and.volume.x.chromosome - mean(dspr.data.reformatted$zhang.radical.score.polarity.and.volume.x.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$zhang.radical.score.polarity.and.volume.x.chromosome,na.rm = TRUE)
dspr.data.reformatted$INDELs.X.chromosome.z.transformed = (dspr.data.reformatted$INDELs.X.chromosome - mean(dspr.data.reformatted$INDELs.X.chromosome,na.rm = TRUE)) / sd(dspr.data.reformatted$INDELs.X.chromosome,na.rm = TRUE)

## Plotting to visualize potential outliers
subset = subset(dspr.data.reformatted, dspr.data.reformatted$sex == "m" & dspr.data.reformatted$mating.regime == "vial")
## Lowest value is X12258 (344) = 28, second lowest is X11287 (168) = 39

outlier.plot=ggplot(subset, aes(x = seq(from = 1, to = 357, by = 1), y = zhang.radical.score.polarity.and.volume.whole.genome)) +
             geom_point(size = 5) +
             theme_bw() +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             axis.text.y = element_text(color="black",size=25, family = "Helvetica"),
             axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank(),
             axis.title.y = element_text(color="black",size=40, family = "Helvetica")) +
             geom_hline(yintercept = mean(subset$zhang.radical.score.polarity.and.volume.whole.genome), size = 2) +
             annotate(geom = "text", x = 344, y = 28, label = "12258", hjust = 1.25, vjust = 1, size = 10) +
             annotate(geom = "text", x = 168, y = 39, label = "11287", hjust = 1.25, vjust = 1, size = 10) +
             ylab("Radical Amino Acid Burden")
ggsave(outlier.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/dspr.outlier.plot.pdf", device = "pdf", units = "in", width = 12, height = 12)



theme_bw()+
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
        legend.position = "none", plot.margin = margin(20, 20, 20, 20)) +
geom_hline(yintercept = 0, size = 1, linetype = "longdash")



## Detecting possible outliers
mean(dspr.data.reformatted$synonymous.variants.whole.genome) # mean = 343.5098
3 * sd(dspr.data.reformatted$synonymous.variants.whole.genome) # 3*sd = 126.7409
lower.bound = mean(dspr.data.reformatted$synonymous.variants.whole.genome) - (3 * sd(dspr.data.reformatted$synonymous.variants.whole.genome))
#No lines have more than 3*SD + mean syn variants but some have less, so lets remove those
dspr.data.reformatted.outliers.excluded = dspr.data.reformatted[dspr.data.reformatted$synonymous.variants.whole.genome > lower.bound,]
dspr.fitness.data.outliers.excluded = dspr.fitness.data[dspr.fitness.data$DSPR.stock.number %in% dspr.data.reformatted.outliers.excluded$DSPR.stock.number,]

# General statistics about deleterious variants
total.variants = (dspr.data.reformatted.outliers.excluded$most.severe.variants.whole.genome +
                 dspr.data.reformatted.outliers.excluded$INDELs.whole.genome +
                 dspr.data.reformatted.outliers.excluded$zhang.radical.score.polarity.and.volume.whole.genome)

print(c(min(total.variants), max(total.variants), mean(total.variants)))

print(c(min(dspr.data.reformatted.outliers.excluded$most.severe.variants.whole.genome),
        max(dspr.data.reformatted.outliers.excluded$most.severe.variants.whole.genome),
        mean(dspr.data.reformatted.outliers.excluded$most.severe.variants.whole.genome)))

print(c(min(dspr.data.reformatted.outliers.excluded$zhang.radical.score.polarity.and.volume.whole.genome),
        max(dspr.data.reformatted.outliers.excluded$zhang.radical.score.polarity.and.volume.whole.genome),
        mean(dspr.data.reformatted.outliers.excluded$zhang.radical.score.polarity.and.volume.whole.genome)))

print(c(min(dspr.data.reformatted.outliers.excluded$INDELs.whole.genome),
        max(dspr.data.reformatted.outliers.excluded$INDELs.whole.genome),
        mean(dspr.data.reformatted.outliers.excluded$INDELs.whole.genome)))



# Merge datasets together
dspr.data.for.model = merge(dspr.data.reformatted.outliers.excluded, dspr.fitness.data.outliers.excluded, by.x = c("DSPR.stock.number", "sex", "mating.regime"), by.y = c("DSPR.stock.number", "sex", "mating.regime"), all = TRUE, sort  = FALSE)

# Update GRM by removing the outlier lines
dspr.GRM = dspr.GRM[rownames(dspr.GRM) %in% dspr.data.reformatted.outliers.excluded$DSPR.stock.number, colnames(dspr.GRM) %in% dspr.data.reformatted.outliers.excluded$DSPR.stock.number]
dspr.GRM.male = dspr.GRM[rownames(dspr.GRM) %in% male.data$DSPR.stock.number, colnames(dspr.GRM) %in% male.data$DSPR.stock.number]
dspr.GRM.female = dspr.GRM[rownames(dspr.GRM) %in% female.data$DSPR.stock.number, colnames(dspr.GRM) %in% female.data$DSPR.stock.number]

dspr.GRM.null = dspr.GRM.null[rownames(dspr.GRM.null) %in% dspr.data.reformatted.outliers.excluded$DSPR.stock.number, colnames(dspr.GRM.null) %in% dspr.data.reformatted.outliers.excluded$DSPR.stock.number]
dspr.GRM.male.null = dspr.GRM.null[rownames(dspr.GRM.null) %in% male.data$DSPR.stock.number, colnames(dspr.GRM.null) %in% male.data$DSPR.stock.number]
dspr.GRM.female.null = dspr.GRM.null[rownames(dspr.GRM.null) %in% female.data$DSPR.stock.number, colnames(dspr.GRM.null) %in% female.data$DSPR.stock.number]

##############################
####     Model fitting    ####
##############################

# Model fitting in lme4 using the line means for fitness
# Full model
dspr.lme4.model = lmer(fitness ~ mating.regime*sex*most.severe.variants.whole.genome + mating.regime*sex*zhang.radical.score.polarity.and.volume.whole.genome + mating.regime*sex*INDELs.whole.genome + (1|DSPR.stock.number), data = dspr.data.reformatted.outliers.excluded)

# Sex specific models
dspr.lme4.model.female = lmer(fitness ~ mating.regime*most.severe.variants.whole.genome+ mating.regime*zhang.radical.score.polarity.and.volume.whole.genome + mating.regime*INDELs.whole.genome + (1|DSPR.stock.number), data = dspr.data.reformatted.outliers.excluded[dspr.data.reformatted.outliers.excluded$sex == "f",])
dspr.lme4.model.male = lmer(fitness ~ mating.regime*most.severe.variants.whole.genome+ mating.regime*zhang.radical.score.polarity.and.volume.whole.genome + mating.regime*INDELs.whole.genome + (1|DSPR.stock.number), data = dspr.data.reformatted.outliers.excluded[dspr.data.reformatted.outliers.excluded$sex == "m",])

# Environment specific models for females
dspr.lme4.model.female.vials = lm(fitness ~ most.severe.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome, data = dspr.data.reformatted.outliers.excluded[dspr.data.reformatted.outliers.excluded$sex == "f" & dspr.data.reformatted.outliers.excluded$mating.regime == "vial",])

dspr.lme4.model.female.vials = lmer(fitness ~ most.severe.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome + (1|DSPR.stock.number), data = dspr.data.reformatted.outliers.excluded[dspr.data.reformatted.outliers.excluded$sex == "f" & dspr.data.reformatted.outliers.excluded$mating.regime == "vial",])
dspr.lme4.model.female.cages = lmer(fitness ~ most.severe.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome + (1|DSPR.stock.number), data = dspr.data.reformatted.outliers.excluded[dspr.data.reformatted.outliers.excluded$sex == "f" & dspr.data.reformatted.outliers.excluded$mating.regime == "cage",])



##  Model fitting -- MCMCglmm models

# This is a code that will run however many runs of the MCMCglmm code and will output all of them to a list
#setCores = 5 # use detectCores() by itself if you want all CPUs

## Getting MCMCglmm to run in parallel
# make the cluster
#cl <- makeCluster(getOption("cl.cores",setCores))

# load the MCMCglmm package within the cluster
#cl.pkg <- clusterEvalQ(cl,library(MCMCglmm))
# import each object that's necessary to run the function

# Setting priors

#prior<-list(R=list(V=diag(0.001,2), nu=1.002), G=list(G1=list(V=diag(0.001, 2), nu=1.002)))

#clusterExport(cl,"prior")
#clusterExport(cl,"female.data")
#clusterExport(cl,"dspr.GRM")
# use parLapply() to execute 5 runs of MCMCglmm(), each with nitt=300000
#DSPR.MCMCglmm.model.female.data.outliers.excluded_5runs<-parLapply(cl=cl,1:15, function(i) {
#  N <- dim(dspr.GRM)[1]
#  i <- rep(1:N,rep(N,N))
#  j <- rep(1:N,N)
#  s <-spMatrix(N,N,i,j,as.vector(dspr.GRM))
#  Ginv<-solve(s)
#  class(Ginv) <- "dgCMatrix"
#  rownames(Ginv) <- Ginv@Dimnames[[1]] <- with(female.data,unique(DSPR.stock.number))
#  model = MCMCglmm(fitness ~ environment*most.severe.variants.whole.genome + environment*zhang.radical.score.polarity.and.volume.whole.genome + environment*INDELs.whole.genome,
#                    random = ~us(environment):DSPR.stock.number,
#                    rcov = ~idh(environment):units,
#                    ginverse=list(DSPR.stock.number=Ginv),
#                    data=female.data,
#                    family = "gaussian",
#                    prior=prior,
#                    nitt = 310000, burnin = 10000, thin = 300)
#                  }
#)
#cluster
#stopCluster(cl)

#saveRDS(DSPR.MCMCglmm.model.female.data.outliers.excluded_5runs, "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/DSPR.MCMCglmm.model.female.data.outliers.excluded_5runs.RDS")

#################
### PLOTTING  ###
#################

## Plotting interactions with ggplot2
require(ggplot2)
require(effects)

dspr.lme4.model = lmer(fitness ~ mating.regime*sex*most.severe.variants.whole.genome + mating.regime*sex*zhang.radical.score.polarity.and.volume.whole.genome + mating.regime*sex*INDELs.whole.genome + (1|DSPR.stock.number), data = dspr.data.reformatted.outliers.excluded)

# Extract interactions into seperate dataframes
# LoF mutations
dspr.effects.most.severe = effect('mating.regime*sex*most.severe.variants.whole.genome', dspr.lme4.model.untransformed,
                          xlevels=list(most.severe.variants.whole.genome = seq(0,12, by = 1)),
                          se=TRUE, confidence.level=.95, typical=mean)

dspr.effects.zhang.score = effect('mating.regime*sex*zhang.radical.score.polarity.and.volume.whole.genome', dspr.lme4.model.untransformed,
                          xlevels=list(zhang.radical.score.polarity.and.volume.whole.genome = seq(190,365, by = 1)),
                          se=TRUE, confidence.level=.95, typical=mean)

dspr.effects.indels = effect('mating.regime*sex*INDELs.whole.genome', dspr.lme4.model.untransformed,
                          xlevels=list(INDELS.whole.genome = seq(75,170, by = 1)),
                          se=TRUE, confidence.level=.95, typical=mean)

#Put data in data frame
dspr.effects.most.severe = as.data.frame(dspr.effects.most.severe)
dspr.effects.zhang.score = as.data.frame(dspr.effects.zhang.score)
dspr.effects.indels = as.data.frame(dspr.effects.indels)

dspr.lof.plot<-ggplot(dspr.effects.most.severe, aes(x=most.severe.variants.whole.genome, y=fit, colour=sex, linetype = mating.regime))+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
              legend.position = "none", plot.margin = margin(20, 20, 20, 20)) +
      geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.001)) +
      scale_fill_manual(values = c("#ff9999", "#9999ff")) +
      geom_line(size=2, aes(color=sex, linetype = mating.regime)) +
      scale_color_manual(values = c("#ff0000", "#0000ff")) +
      scale_x_continuous(breaks=seq(0,12,3)) +
      coord_cartesian(xlim=c(0,12)) +
      geom_hline(yintercept = 0, size = 1, linetype = "longdash")

dspr.zhang.plot<-ggplot(dspr.effects.zhang.score, aes(x=zhang.radical.score.polarity.and.volume.whole.genome, y=fit, colour=sex, linetype = mating.regime))+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
              legend.position = "none", plot.margin = margin(20, 20, 20, 20)) +
      geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.001)) +
      scale_fill_manual(values = c("#ff9999", "#9999ff")) +
      geom_line(size=2, aes(color=sex, linetype = mating.regime)) +
      scale_color_manual(values = c("#ff0000", "#0000ff")) +
      scale_x_continuous(breaks=seq(190,370,60)) +
      geom_hline(yintercept = 0, size = 1, linetype = "longdash")

dspr.indels.plot<-ggplot(dspr.effects.indels, aes(x=INDELs.whole.genome, y=fit, colour=sex, linetype = mating.regime))+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
              legend.position = "none", plot.margin = margin(20, 20, 20, 20)) +
      geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.001)) +
      scale_fill_manual(values = c("#ff9999", "#9999ff")) +
      geom_line(size=3, aes(color=sex, linetype = mating.regime)) +
      scale_color_manual(values = c("#ff0000", "#0000ff")) +
      scale_x_continuous(breaks=seq(75,175,25)) +
      scale_y_continuous( limits = c(-0.1,0.1), expand = c(0,0)) +
      coord_cartesian(xlim=c(75,175)) +
      geom_hline(yintercept = 0, size = 1, linetype = "longdash")

ggsave(dspr.lof.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/dspr.lof.plot.pdf", device = "pdf", units = "in", width = 12, height = 24)
ggsave(dspr.zhang.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/dspr.zhang.plot.pdf", device = "pdf", units = "in", width = 12, height = 24)
ggsave(dspr.indels.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/dspr.indels.plot.pdf", device = "pdf", units = "in", width = 12, height = 24)


# Plotting seperately for each mating regime


theme_bw()+
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
        legend.position = "none") +
#scale_fill_grey() +
geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.25)) +
scale_fill_manual(values = c("#ff9999", "#9999ff")) +
geom_line(size=3, aes(color=sex)) +
scale_color_manual(values = c("#ff0000", "#0000ff")) +
scale_x_continuous(breaks=seq(100,170,10)) +
geom_hline(yintercept = 0, size = 1, linetype = "longdash")




### Full model
## Make a dataframe of effects for each term in the model
dspr.data.model.effects.most.severe = as.data.frame(effect("mating.regime*sex*most.severe.variants.whole.genome",dspr.lme4.model.untransformed))
dspr.data.model.effects.zhang.index = as.data.frame(effect("mating.regime*sex*zhang.radical.score.polarity.and.volume.whole.genome",dspr.lme4.model.untransformed))
dspr.data.model.effects.indels = as.data.frame(effect("mating.regime*sex*INDELs.whole.genome",dspr.lme4.model.untransformed))
names(dspr.data.model.effects.most.severe)[3] = "x.axis.scale"
names(dspr.data.model.effects.zhang.index)[3] = "x.axis.scale"
names(dspr.data.model.effects.indels)[3] = "x.axis.scale"
dspr.data.model.effects.most.severe$trait = "most.severe.variants.whole.genome"
dspr.data.model.effects.zhang.index$trait = "zhang.radical.score.polarity.and.volume.whole.genome"
dspr.data.model.effects.indels$trait = "INDELS.whole.genome"
# Join all data into a single data frame
dspr.data.model.effects = rbind(dspr.data.model.effects.most.severe,
                                dspr.data.model.effects.zhang.index,
                                dspr.data.model.effects.indels)

dspr.data.model.effects.cage = dspr.data.model.effects[dspr.data.model.effects$mating.regime == "cage", ]
dspr.data.model.effects.vial = dspr.data.model.effects[dspr.data.model.effects$mating.regime == "vial", ]

theme_set(theme_bw())




plot_model(dspr.data.model.effects.cage, type = "pred", terms = c("sex", "trait"))
plot_model(dspr.lme4.model.untransformed, type = "int")

pdf(file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/modeloutput.pdf", width = 12, height = 24)
dspr.model.plot.cage = ggplot(dspr.data.model.effects, aes(y = fit,x = x.axis.scale, colour = sex, linetype = mating.regime)) +
    geom_line() +
    ## colour=NA suppresses edges of the ribbon
    geom_ribbon(colour=NA,alpha=0.1, aes(ymin=lower,ymax=upper)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~trait, scales = "free_x")
dspr.model.plot.cage

dspr.model.plot.vial = ggplot(as.data.frame(dspr.data.model.effects.vial),aes(y = fit,x = x.axis.scale, colour = sex)) +
    geom_line() +
    ## colour=NA suppresses edges of the ribbon
    geom_ribbon(colour=NA,alpha=0.1, aes(ymin=lower,ymax=upper)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~trait, scales = "free_x")
dspr.model.plot.vial

theme_set(theme_bw())
dspr.model.plot = ggplot(as.data.frame(dspr.data.model.effects),aes(y = fit,x = x.axis.scale, colour = sex, fill = environment)) +
    geom_line(aes(linetype = environment)) +
    ## colour=NA suppresses edges of the ribbon
    geom_ribbon(colour=NA,alpha=0.1, aes(ymin=lower,ymax=upper)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~trait, scales = "free_x")
dspr.model.plot



###############################################
######### Analysis with Ruzicka Data  #########
###############################################

# Loading in phenotypic data
ruzicka.data = read.delim("/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Ruzicka.Data/pheno.data.with.variants.Feb24.2020.txt", header = TRUE)
# Loading in genetic relatedness matrix matrix
lhm.GRM = as.matrix(read.delim("/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Ruzicka.Data/lhm.kinship.matrix/lhm.kinship.matrix.Feb26.txt", head = TRUE))
#Preparing GRM
# Make the GRM symmetric
lowerTriangle(lhm.GRM) <- upperTriangle(lhm.GRM, byrow=TRUE) #This is a function that will convert the upper triangle to the lower triangle
# Normalize the GRM to
# lhm.GRM.smoothed = as.matrix(cor.smooth(lhm.GRM)) #This is a function that will make the GRM positive definite
lhm.GRM = make.positive.definite(lhm.GRM) #This is a function that will make the GRM positive definite



## Formatting phenotype data for model
ruzicka.data.male.tmp = ruzicka.data[, -c(2,4:10)]
ruzicka.data.female.tmp = ruzicka.data[, -c(2,3,5:10)]
names(ruzicka.data.male.tmp)[2] = "fitness"
names(ruzicka.data.female.tmp)[2] = "fitness"
ruzicka.data.male.tmp$sex = "m"
ruzicka.data.female.tmp$sex = "f"
ruzicka.data.reformatted = rbind(ruzicka.data.male.tmp, ruzicka.data.female.tmp)

##  Remove lines that did not have genotypic information
row.has.na = apply(ruzicka.data.reformatted, 1, function(x) {any(is.na(x))})
ruzicka.data.reformatted.filtered = ruzicka.data.reformatted[!row.has.na,]

## Scale each trait by performing a Z transformation to each
ruzicka.data.reformatted.filtered$most.severe.variants.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$most.severe.variants.whole.genome - mean(ruzicka.data.reformatted.filtered$most.severe.variants.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$most.severe.variants.whole.genome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome - mean(ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome, na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome, na.rm=TRUE)
ruzicka.data.reformatted.filtered$missense.variants.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$missense.variants.whole.genome - mean(ruzicka.data.reformatted.filtered$missense.variants.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$missense.variants.whole.genome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$stop.gained.variants.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$stop.gained.variants.whole.genome - mean(ruzicka.data.reformatted.filtered$stop.gained.variants.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$stop.gained.variants.whole.genome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$stop.lost.variants.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$stop.lost.variants.whole.genome - mean(ruzicka.data.reformatted.filtered$stop.lost.variants.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$stop.lost.variants.whole.genome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$splice.disrupting.variants.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$splice.disrupting.variants.whole.genome - mean(ruzicka.data.reformatted.filtered$splice.disrupting.variants.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$splice.disrupting.variants.whole.genome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.whole.genome - mean(ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.whole.genome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$INDELS.whole.genome.z.transformed = (ruzicka.data.reformatted.filtered$INDELS.whole.genome - mean(ruzicka.data.reformatted.filtered$INDELS.whole.genome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$INDELS.whole.genome,na.rm = TRUE)

ruzicka.data.reformatted.filtered$most.severe.variants.autosome.z.transformed = (ruzicka.data.reformatted.filtered$most.severe.variants.autosome - mean(ruzicka.data.reformatted.filtered$most.severe.variants.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$most.severe.variants.autosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$synonymous.variants.autosome.z.transformed = (ruzicka.data.reformatted.filtered$synonymous.variants.autosome - mean(ruzicka.data.reformatted.filtered$synonymous.variants.autosome, na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$synonymous.variants.autosome, na.rm=TRUE)
ruzicka.data.reformatted.filtered$missense.variants.autosome.z.transformed = (ruzicka.data.reformatted.filtered$missense.variants.autosome - mean(ruzicka.data.reformatted.filtered$missense.variants.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$missense.variants.autosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$stop.gained.variants.autosome.z.transformed = (ruzicka.data.reformatted.filtered$stop.gained.variants.autosome - mean(ruzicka.data.reformatted.filtered$stop.gained.variants.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$stop.gained.variants.autosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$stop.lost.variants.autosome.z.transformed = (ruzicka.data.reformatted.filtered$stop.lost.variants.autosome - mean(ruzicka.data.reformatted.filtered$stop.lost.variants.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$stop.lost.variants.autosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$splice.disrupting.variants.autosome.z.transformed = (ruzicka.data.reformatted.filtered$splice.disrupting.variants.autosome - mean(ruzicka.data.reformatted.filtered$splice.disrupting.variants.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$splice.disrupting.variants.autosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.autosome.z.transformed = (ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.autosome - mean(ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.autosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$INDELS.autosome.z.transformed = (ruzicka.data.reformatted.filtered$INDELS.autosome - mean(ruzicka.data.reformatted.filtered$INDELS.autosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$INDELS.autosome,na.rm = TRUE)

ruzicka.data.reformatted.filtered$most.severe.variants.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$most.severe.variants.X.chromosome - mean(ruzicka.data.reformatted.filtered$most.severe.variants.X.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$most.severe.variants.X.chromosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$synonymous.variants.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$synonymous.variants.X.chromosome - mean(ruzicka.data.reformatted.filtered$synonymous.variants.X.chromosome, na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$synonymous.variants.X.chromosome, na.rm=TRUE)
ruzicka.data.reformatted.filtered$missense.variants.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$missense.variants.X.chromosome - mean(ruzicka.data.reformatted.filtered$missense.variants.X.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$missense.variants.X.chromosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$stop.gained.variants.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$stop.gained.variants.X.chromosome - mean(ruzicka.data.reformatted.filtered$stop.gained.variants.X.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$stop.gained.variants.X.chromosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$stop.lost.variants.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$stop.lost.variants.X.chromosome - mean(ruzicka.data.reformatted.filtered$stop.lost.variants.X.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$stop.lost.variants.X.chromosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$splice.disrupting.variants.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$splice.disrupting.variants.X.chromosome - mean(ruzicka.data.reformatted.filtered$splice.disrupting.variants.X.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$splice.disrupting.variants.X.chromosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.x.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.x.chromosome - mean(ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.x.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$zhang.radical.score.polarity.and.volume.x.chromosome,na.rm = TRUE)
ruzicka.data.reformatted.filtered$INDELS.X.chromosome.z.transformed = (ruzicka.data.reformatted.filtered$INDELS.X.chromosome - mean(ruzicka.data.reformatted.filtered$INDELS.X.chromosome,na.rm = TRUE)) / sd(ruzicka.data.reformatted.filtered$INDELS.X.chromosome,na.rm = TRUE)


# Based on visual inspection, lines H072 and H093 seem like outliers
# H072 and H093 are the only lines that are greater than 3 SDs away from the mean with respect to synonymous variants (and others)

# A more formal assessment of outliers is a Grubb's test which tests a distribution of values for a evidence that a single
# point is an outlier. The advice that I found was that I should run the analysis, remove the offending point, rerun until there is no evidence of an oitlier
#grubbs.test(ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome)
# Remove largest value as result of this test
#ruzicka.data.reformatted.filtered = ruzicka.data.reformatted.filtered[ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome < 62371,]
#grubbs.test(ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome)
#ruzicka.data.reformatted.filtered = ruzicka.data.reformatted.filtered[ruzicka.data.reformatted.filtered$synonymous.variants.whole.genome < 60400,]

#Plotting outliers for supplemental
subset = subset(ruzicka.data.reformatted.filtered, ruzicka.data.reformatted.filtered$sex == "m")
## Highest value is H093 ( 8617 ) = 89, second highest is H072 (8070) = 70

outlier.plot=ggplot(subset, aes(x = seq(from = 1, to = 220, by = 1), y = zhang.radical.score.polarity.and.volume.whole.genome)) +
                    geom_point(size = 5) +
                    theme_bw() +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    axis.text.y = element_text(color="black",size=25, family = "Helvetica"),
                    axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank(),
                    axis.title.y = element_text(color="black",size=40, family = "Helvetica")) +
                    geom_hline(yintercept = mean(subset$zhang.radical.score.polarity.and.volume.whole.genome), size = 2) +
                    annotate(geom = "text", x = 89, y = 9649, label = "H093", hjust = 1.5, vjust = 1, size = 10) +
                    annotate(geom = "text", x = 70, y = 9009, label = "H072", hjust = 1.5, vjust = 1, size = 10) +
                    ylab("Radical Amino Acid Burden")
#outlier.plot
ggsave(outlier.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/ruzicka.outlier.plot.pdf", device = "pdf", units = "in", width = 12, height = 12)

# Exclude lines H072 and H093 as these are the two lines that are likely outliers
ruzicka.data.reformatted.filtered = ruzicka.data.reformatted.filtered[ruzicka.data.reformatted.filtered$IID != "H072" & ruzicka.data.reformatted.filtered$IID != "H093", ]
# Update GRM by removing the outlier lines
lhm.GRM = lhm.GRM[rownames(lhm.GRM) %in% ruzicka.data.reformatted.filtered$IID, colnames(lhm.GRM) %in% ruzicka.data.reformatted.filtered$IID]



##############################
####     Model fitting    ####
##############################

### Testing the effect of the Kinship matrix ###
# Make an identiy matrix
identity.matrix = diag(1, nrow = nrow(lhm.GRM), ncol =ncol(lhm.GRM))
rownames(identity.matrix) = rownames(lhm.GRM)
colnames(identity.matrix) = colnames(lhm.GRM)
# Fit a kinship matrix and identiy matrix as the (co)variance matrix in seperate models
ruzicka.lmekin.model.kinship.matrix = lmekin(fitness ~ (1|IID), data = ruzicka.data.reformatted.filtered, varlist = list(lhm.GRM), method = "ML")
ruzicka.lmekin.model.identity.matrix = lmekin(fitness ~ (1|IID), data = ruzicka.data.reformatted.filtered, varlist = list(identity.matrix), method = "ML")

# Model fitting in lme4
ruzicka.lme4.model = lmer(fitness ~ sex*most.severe.variants.whole.genome + sex*zhang.radical.score.polarity.and.volume.whole.genome + sex*INDELS.whole.genome + (1|IID), data = ruzicka.data.reformatted.filtered)


## Plotting interactions
library(effects)
# Extract interactions into seperate dataframes
# LoF mutations
ruzicka.effects.most.severe = effect('sex*most.severe.variants.whole.genome', ruzicka.lme4.model,
                          xlevels=list(most.severe.variants.whole.genome = seq(110,160, by = 1)),
                          se=TRUE, confidence.level=.95, typical=mean)

ruzicka.effects.zhang.score = effect('sex*zhang.radical.score.polarity.and.volume.whole.genome', ruzicka.lme4.model,
                          xlevels=list(zhang.radical.score.polarity.and.volume.whole.genome = seq(8300,8700, by = 1)),
                          se=TRUE, confidence.level=.95, typical=mean)

ruzicka.effects.indels = effect('sex*INDELS.whole.genome', ruzicka.lme4.model,
                          xlevels=list(INDELS.whole.genome = seq(180,320, by = 1)),
                          se=TRUE, confidence.level=.95, typical=mean)

#Put data in data frame
ruzicka.effects.most.severe<-as.data.frame(ruzicka.effects.most.severe)
ruzicka.effects.zhang.score<-as.data.frame(ruzicka.effects.zhang.score)
ruzicka.effects.indels<-as.data.frame(ruzicka.effects.indels)

#Create factors of the interaction variables
#ruzicka.effects.most.severe$most.severe = factor(ruzicka.effects.most.severe$most.severe.whole.genome,levels=seq(110,160,1))
#ruzicka.effects.zhang.score$zhang.radical.score.polarity.and.volume.whole.genome = numeric(ruzicka.effects.zhang.score$zhang.radical.score.polarity.and.volume.whole.genome,levels=seq(8300,8700,1))
#ruzicka.effects.indels$INDELS.whole.genome = factor(ruzicka.effects.indels$INDELS.whole.genome, levels=seq(180,320,1))

ruzicka.lof.plot<-ggplot(ruzicka.effects.most.severe, aes(x=most.severe.variants.whole.genome, y=fit, group=sex))+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
              legend.position = "none") +
      #scale_fill_grey() +
      geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.25)) +
      scale_fill_manual(values = c("#ff9999", "#9999ff")) +
      geom_line(size=3, aes(color=sex)) +
      scale_color_manual(values = c("#ff0000", "#0000ff")) +
      scale_x_continuous(breaks=seq(100,170,10)) +
      geom_hline(yintercept = 0, size = 1, linetype = "longdash")


ruzicka.zhang.plot<-ggplot(ruzicka.effects.zhang.score, aes(x=zhang.radical.score.polarity.and.volume.whole.genome, y=fit, group=sex))+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
              legend.position = "none") +
      geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.25)) +
      scale_fill_manual(values = c("#ff9999", "#9999ff")) +
      geom_line(size=2, aes(color=sex)) +
      scale_color_manual(values = c("#ff0000", "#0000ff")) +
      scale_x_continuous(breaks=seq(8000,9000,100)) +
      geom_hline(yintercept = 0, size = 1, linetype = "longdash")

ruzicka.indels.plot<-ggplot(ruzicka.effects.indels, aes(x=INDELS.whole.genome, y=fit, group=sex))+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.text = element_text(face="bold", color="black",size=40, family = "Helvetica"), axis.title=element_blank(),
              legend.position = "none") +
      geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = sex, color = sex, alpha = 0.25)) +
      scale_fill_manual(values = c("#ff9999", "#9999ff")) +
      geom_line(size=2, aes(color=sex)) +
      scale_color_manual(values = c("#ff0000", "#0000ff")) +
      scale_x_continuous(breaks=seq(180,320,20)) +
      geom_hline(yintercept = 0, size = 1, linetype = "longdash")

# Saving plots for Ruzicka DATA
ggsave(ruzicka.lof.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/ruzicka.lof.plot.pdf", device = "pdf", units = "in", width = 12, height = 24)
ggsave(ruzicka.zhang.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/ruzicka.zhang.plot.pdf", device = "pdf", units = "in", width = 12, height = 24)
ggsave(ruzicka.indels.plot, file = "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/Plots/ruzicka.indels.plot.pdf", device = "pdf", units = "in", width = 12, height = 24)







# Plotting Ruzicka data
effects.most.severe = effect("sex*most.severe.variants.whole.genome",ruzicka.lme4.model.untransformed , KR = T)
plot(effects.most.severe.test)
effects.radical.score = effect("sex*zhang.radical.score.polarity.and.volume.whole.genome",ruzicka.lme4.model.untransformed , KR = T)
plot(effects.radical.score)
effects.indels = effect("sex*INDELS.whole.genome",ruzicka.lme4.model.untransformed , KR = T)
plot(effects.indels)

##########
## Model fitting with MCMCglmm  -- This ended up not being a route we went down
# This is a code that will run however many runs of the MCMCglmm code and will output all of them to a list

#setCores = 5 # use detectCores() by itself if you want all CPUs

## Getting MCMCglmm to run in parallel
# make the cluster
#cl <- makeCluster(getOption("cl.cores",setCores))

# load the MCMCglmm package within the cluster
#cl.pkg <- clusterEvalQ(cl,library(MCMCglmm))
# import each object that's necessary to run the function

#prior<-list(R=list(V=diag(2), nu=1.002), G=list(G1=list(V=diag(2), nu=1.002)))


#clusterExport(cl,"prior")
#clusterExport(cl,"ruzicka.data.reformatted.filtered")
#clusterExport(cl,"lhm.GRM")
# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
#ruzicka.MCMCglmm.full.model_5runs<-parLapply(cl=cl,1:10, function(i) {
#N <- dim(lhm.GRM)[1]
#i <- rep(1:N,rep(N,N))
#j <- rep(1:N,N)
#s <-spMatrix(N,N,i,j,as.vector(lhm.GRM))
#Ginv<-solve(s)
#class(Ginv) <- "dgCMatrix"
#rownames(Ginv) <- Ginv@Dimnames[[1]] <- with(ruzicka.data.reformatted.filtered,unique(IID))

# Setting priors
#model = MCMCglmm(fitness ~ sex*most.severe.variants.whole.genome.z.transformed + sex*zhang.radical.score.polarity.and.volume.whole.genome.z.transformed + sex*INDELS.whole.genome.z.transformed,
#                random = ~us(sex):IID,
#                rcov = ~idh(sex):units,
#                ginverse=list(IID=Ginv),
#                data=ruzicka.data.reformatted.filtered,
#                family = "gaussian",
#                prior=prior,
#                nitt = 301000, burnin = 10000, thin = 300)
#})

#cluster
#stopCluster(cl)
#saveRDS(ruzicka.MCMCglmm.full.model_5runs, "/plas1/amardeep.singh/GenomicsMaleFemaleArchitecture/Mutation.Load.Project/ruzicka.MCMCglmm.full.model_5runs.RDS")


#################
### PLOTTING  ###
#################

## Make a dataframe of effects for each term in the model
ruzicka.data.model.effects.most.severe = as.data.frame(effect(c("sex*most.severe.variants.whole.genome"),ruzicka.lme4.model.untransformed))
ruzicka.data.model.effects.zhang.index = as.data.frame(effect(c("sex*zhang.radical.score.polarity.and.volume.whole.genome"),ruzicka.lme4.model.untransformed))
ruzicka.data.model.effects.indels = as.data.frame(effect(c("sex*INDELS.whole.genome"),ruzicka.lme4.model.untransformed))
names(ruzicka.data.model.effects.most.severe)[2] = "x.axis.scale"
names(ruzicka.data.model.effects.zhang.index)[2] = "x.axis.scale"
names(ruzicka.data.model.effects.indels)[2] = "x.axis.scale"
ruzicka.data.model.effects.most.severe$trait = "most.severe.variants.whole.genome"
ruzicka.data.model.effects.zhang.index$trait = "zhang.radical.score.polarity.and.volume.whole.genome"
ruzicka.data.model.effects.indels$trait = "INDELS.whole.genome"
# Join all data into a single data frame
ruzicka.data.model.effects = rbind(ruzicka.data.model.effects.most.severe,ruzicka.data.model.effects.zhang.index,ruzicka.data.model.effects.indels)

theme_set(theme_bw())
ruzicka.model.plot = ggplot(as.data.frame(ruzicka.data.model.effects),aes(y = fit,x = x.axis.scale,colour=sex,fill=sex)) +
    geom_line() +
    ## colour=NA suppresses edges of the ribbon
    geom_ribbon(colour=NA,alpha=0.1, aes(ymin=lower,ymax=upper)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~trait, scales = "free_x")
ruzicka.model.plot


########################################################################################################################
########################################################################################################################









## Testing Grounds


# Female Models
dspr.lme4.model.untransformed.female
dspr.data.model.effects.most.severe = effect("environment*most.severe.variants.whole.genome",dspr.lme4.model.untransformed.female )
dspr.data.model.effects.zhang.index = effect("environment*zhang.radical.score.polarity.and.volume.whole.genome",dspr.lme4.model.untransformed.female)
dspr.data.model.effects.indels = effect("environment*INDELs.whole.genome",dspr.lme4.model.untransformed.female)








# Full model



names(dspr.data.model.effects.most.severe)[3] = "x.axis.scale"
names(dspr.data.model.effects.zhang.index)[3] = "x.axis.scale"
names(dspr.data.model.effects.indels)[3] = "x.axis.scale"

dspr.data.model.effects.most.severe$trait = "most.severe.variants.whole.genome"
dspr.data.model.effects.zhang.index$trait = "zhang.radical.score.polarity.and.volume.whole.genome"
dspr.data.model.effects.indels$trait = "INDELs.whole.genome"


dspr.data.model.effects = rbind(dspr.data.model.effects.most.severe,dspr.data.model.effects.zhang.index,dspr.data.model.effects.indels)

theme_set(theme_bw())

plot = ggplot(as.data.frame(dspr.data.model.effects),aes(y = fit,x = x.axis.scale,colour=sex,fill=sex, linetype = environment)) +
    geom_line() +
    ## colour=NA suppresses edges of the ribbon
    geom_ribbon(colour=NA,alpha=0.1, aes(ymin=lower,ymax=upper)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~trait)



# Female model
dspr.data.model.effects.most.severe = as.data.frame(effect(c("environment*most.severe.variants.whole.genome.z.transformed"),dspr.data.model.female.environment ))
dspr.data.model.effects.zhang.index = as.data.frame(effect(c("environment*zhang.radical.score.polarity.and.volume.whole.genome.z.transformed"),dspr.data.model.female.environment ))
dspr.data.model.effects.indels = as.data.frame(effect(c("environment*INDELs.whole.genome.z.transformed"),dspr.data.model.female.environment ))

names(dspr.data.model.effects.most.severe)[2] = "x.axis.scale"
names(dspr.data.model.effects.zhang.index)[2] = "x.axis.scale"
names(dspr.data.model.effects.indels)[2] = "x.axis.scale"

dspr.data.model.effects.most.severe$trait = "most.severe.variants.whole.genome"
dspr.data.model.effects.zhang.index$trait = "zhang.radical.score.polarity.and.volume.whole.genome"
dspr.data.model.effects.indels$trait = "INDELs.whole.genome"


dspr.data.model.effects = rbind(dspr.data.model.effects.most.severe,dspr.data.model.effects.zhang.index,dspr.data.model.effects.indels)

    theme_set(theme_bw())

plot = ggplot(as.data.frame(dspr.data.model.effects),aes(y = fit,x = x.axis.scale,colour=environment,fill=environment)) +
    geom_line() +
    ## colour=NA suppresses edges of the ribbon
    geom_ribbon(colour=NA,alpha=0.1, aes(ymin=lower,ymax=upper)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~trait)







# Model fitting
# Assessing the significance of the GRM as a random effect
dspr.lmekin.relatedness.only = lmekin(fitness ~ 1 + (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM), method = "REML")
dspr.lmekin.relatedness.only.null = lmekin(fitness ~ 1 + (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM.null), method = "REML")

lmekin.anova(dspr.lmekin.relatedness.only,dspr.lmekin.relatedness.only.null)

# NOTE: I don't know how to compare a model with and without a random effect using nlme. However, Ben Bolker suggested using AIC to compare a model fit with lme4 (I'm using nlme) fit with ML
#       and NOT REML can be compared to an LM model

AICc(dspr.lmekin.relatedness.only, dspr.lmekin.relatedness.only.null)

#dspr.lmekin.model.sex.and.environment = lmekin(fitness ~ sex + environment + most.severe.variants.whole.genome + missense.variants.whole.genome + stop.gained.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome +
#    sex*most.severe.variants.whole.genome + sex*missense.variants.whole.genome + sex*stop.gained.variants.whole.genome + sex*zhang.radical.score.polarity.and.volume.whole.genome + sex*INDELs.whole.genome +
#    environment*most.severe.variants.whole.genome + environment*missense.variants.whole.genome + environment*stop.gained.variants.whole.genome + environment*zhang.radical.score.polarity.and.volume.whole.genome + environment*INDELs.whole.genome +
#    sex*environment*most.severe.variants.whole.genome + sex*environment*missense.variants.whole.genome + sex*environment*stop.gained.variants.whole.genome + sex*environment*zhang.radical.score.polarity.and.volume.whole.genome + sex*environment*INDELs.whole.genome +
#    (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM), method = "REML")

#dspr.lmekin.model.sex = lmekin(fitness ~ sex + most.severe.variants.whole.genome + missense.variants.whole.genome + stop.gained.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome +
#    sex*most.severe.variants.whole.genome + sex*missense.variants.whole.genome + sex*stop.gained.variants.whole.genome + sex*zhang.radical.score.polarity.and.volume.whole.genome + sex*INDELs.whole.genome +
#    (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM), method = "REML")

#dspr.lmekin.model.environment = lmekin(fitness ~ environment + most.severe.variants.whole.genome + missense.variants.whole.genome + stop.gained.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome +
#    environment*most.severe.variants.whole.genome + environment*missense.variants.whole.genome + environment*stop.gained.variants.whole.genome + environment*zhang.radical.score.polarity.and.volume.whole.genome + environment*INDELs.whole.genome +
#    (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM), method = "REML")

#dspr.lmekin.model = lmekin(fitness ~ most.severe.variants.whole.genome + missense.variants.whole.genome + stop.gained.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome +
#      (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM), method = "REML")

#dspr.lmekin.model.drop.indels = lmekin(fitness ~ most.severe.variants.whole.genome + missense.variants.whole.genome + stop.gained.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome  +
#        (1|DSPR.stock.number), data = dspr.data.reformatted, varlist = list(dspr.GRM), method = "REML")



#AICc(dspr.lmekin.model.sex.and.environment, dspr.lmekin.model.sex)







##
