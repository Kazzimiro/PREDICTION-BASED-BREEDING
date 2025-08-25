###################################################
#HOMEWORK - HS 541 Plant Breeding Methods - 6TH QUESTION
#  DISCENTE: JOSÉ ARTUR DE OLIVEIRA CASIMIRO
#  DOCENTES: ROBERTO FRITSCHE-NETO E JULIO CESAR DoVALE

#contato: artur.casimiro@alu.ufc.br

#atualização: 24/08/2025

###################################################
# loading packages
require(foreach)
require(doParallel)
require(doMC)
library(car)
library(ggplot2)
library(breedR)

# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1) # type the number of cores you want to use
getDoParWorkers()

#loading phenotypes
pheno <- readRDS("pheno")
head(pheno)
dim(pheno)
str(pheno)

# reorganizing the file
pheno <- pheno[,c(1:9, 13:14, 10:12)]
head(pheno)

# index for traits
traits <- colnames(pheno)[12:14]
  
# estimate de phenotypic correlation
pheno.cor <- round(cor(pheno[,traits], use = "pairwise.complete.obs", method = "pearson"), 2)
corrplot::corrplot(pheno.cor, method = 'number', type = "lower", diag = F, col = c("red", "orange", "green", "blue"))

# reorganize de data
pheno.melted <- reshape2::melt(pheno, measure.vars = traits)
head(pheno.melted)  
tail(pheno.melted)

# load the GRM
Ga <- readRDS("Ga")

# running all GS single-traits in parallel 
results.st <- foreach(i = 1:length(traits), 
                        .packages = c("breedR", "car"), 
                        .combine = "rbind",
                        .export = c("remlf90", "outlierTest"),
                        .multicombine = TRUE, 
                        .errorhandling = "remove",
                        .verbose = TRUE    
  ) %do% {

    # subset the data  
    smpl <- droplevels.data.frame(pheno.melted[pheno.melted$variable == traits[i],])
    
    # outlier detection and elimination
    fit <- lm(value ~ rep + gid + N + gid:N, data = smpl)
    outlier <- names(outlierTest(fit)$p)
    smpl[outlier, "value"] <- NA
    
    Za <- model.matrix(~ -1 + gid, data = smpl)
    colnames(Za) <- gsub("gid", "", colnames(Za), fixed = T)

    # MME using only the classical experimental design, and gid as random
    sol <- remlf90(fixed = value ~ rep + N, 
                   generic = list(BV = list(Za, Ga)), # + GN 
                   data = smpl, 
                   method = "em")
    
    # reorganizing BLUPS
    BLUPS <- data.frame(gid = colnames(Za),
                        trait = traits[i],
                        GEBV = sol$ranef$BV[[1]][,1])

}

head(results.st)
tail(results.st)

# GEBVs per trait
GEBVs <- tapply(results.st$GEBV, list(results.st$gid, results.st$trait), mean)
head(GEBVs)

######################## single step (ss) multi-trait (MT) GS ################

# covariance matrices - example with 2 traits
covg <- matrix(c(1, cov(pheno[,12], pheno[,13], use = "pairwise.complete.obs"), 
                 cov(pheno[,12], pheno[,13], use = "pairwise.complete.obs"), 1), 2, 2)
colnames(covg) <- rownames(covg) <- colnames(pheno)[12:13]
covg

covr <- matrix(c(1, 0, 0, 1), 2, 2)
colnames(covr) <- rownames(covr) <- colnames(pheno)[12:13]
covr

sol <- remlf90(fixed = cbind(pheno[,12], pheno[,13]) ~ N + rep, 
               generic = list(GEBV = list(Za, Ga, var.ini = covg)), 
               data = pheno, 
               method = "em", 
               var.ini = list(covr))

# G and R varcomp via MT-GBLUP
sol$var

## predicted values for the test set
blupsA <- as.matrix(sol$ranef[[1]][[1]])
rownames(blupsA) <- colnames(Za)
blupsB <- as.matrix(sol$ranef[[1]][[2]])
rownames(blupsB) <- colnames(Za)

# heritabilities
(hm.A <- sol$var[[1]][1,1]/sum(sol$var[[1]][1,1], sol$var[[2]][1,1]))
(hm.B <- sol$var[[1]][2,2]/sum(sol$var[[1]][2,2], sol$var[[2]][2,2]))
(LL <- logLik(sol)[1])
###observa-se baixa herdabilidade entre os fatores

GEBV.MT = data.frame(
  gid = colnames(Za),
  GEBV_A = blupsA[,1], 
  GEBV_B = blupsB[,1])

# the new correlation between the traits
cor(GEBV.MT[,2:3])
# correlation among traits via MT-GBLUP
cov2cor(sol$var$GEBV)
# the residual correlation
cov2cor(sol$var$Residual)

# So, the only thing now is just weight the GEBV and obtain the SI
# define the economic weights per trait
ecoW <- c(1, 1) # in this case two times more for grain yield
# then, the vector of SI per genotype
MT <- as.matrix(GEBV.MT[,2:3]) %*% ecoW  
head(MT)

########################### Selection indices ##################################
# phenotypic covariance between traits
P <- cov(pheno[, 12:13], use = "pairwise", method = "pearson")

# genetic covariance between traits
G <- cov(GEBVs[,2:3])

# Smith-Hazel
# define the economic weights per trait
ecoW <- c(1, 1) # in this case two times more for grain yield

# then, the selection weights per trait
(b.sh <- solve(as.matrix(P)) %*% as.matrix(G) %*% as.matrix(ecoW))

#[1,1]#
##SDM 0.0001224927
##SRA 0.0165739678##

#[1,2]#
##SDM 0.0001380601##
##SRA 0.0197787274##

#
# Finally, the vector of SI per genotype
SH <- GEBVs[,2:3] %*% b.sh  
head(SH)

# Pasek-Baker
# define the desired genetic gains per traits in genetic standard deviations
desired <- c(1, 1)   .
# G correlation matrix x desired genetic gains in standard deviations
(b.pb <- MASS::ginv(cov(scale(GEBVs[,2:3]))) %*% as.matrix(desired))


#   [1,1]  #
##[1,] 0.5331413##
##[2,] 0.5331413##

#  [1,2]   #
##[1,] -3.222013##
##[2,]  4.821437##


#  [1,3]   #
##[1,] -6.977167##
##[2,]  9.109733##


PB <- as.matrix(scale(GEBVs[,2:3])) %*% b.pb
head(PB)

# correlation between GEBVs and selection indices
GEBVs <- as.data.frame(GEBVs[,2:3])
GEBVs$SH <- SH
GEBVs$PB <- PB
GEBVs$MT <- MT
head(GEBVs)
cor(GEBVs)



##for ecow (1, 1)

#          PB          MT            SH
#1    0.4248290  0.00937835  0.0002318921
#10   3.2073344  0.07455769 -0.0008982987
#11  -1.3724195 -0.00361111 -0.0008058520
#110  0.5833876 -0.01089990 -0.0003332033
#112  1.1910198  0.02431674  0.0002729751
#113  2.3345741  0.04627886  0.0001196571



##for ecow (1, 2)

#             SH         PB          MT
# 1    0.0002318921  0.8784754  0.01039056
# 10  -0.0008982987  7.4579069  0.07073444
# 11  -0.0008058520 -6.3008132 -0.00827387
# 110 -0.0003332033  4.1682038 -0.00638074
# 112  0.0002729751  1.8903961  0.03255225
# 113  0.0001196571  5.5823891  0.05050219



# A principal conclusão é que a alteração dos pesos muda a prioridade da seleção.

# Na análise de linha de base (com pesos iguais de 1,1), o índice MT ranqueou os genótipos de uma maneira.

# Ao dobrar o peso da característica SRA (pesos de 1,2), a pontuação do índice MT para cada genótipo se altera.

# Isso acontece porque o ranqueamento é agora um reflexo direto da nova prioridade de seleção.

# A alteração dos pesos econômicos influencia a lista final de genótipos selecionados.










# let's see which method provides the best KPIs
# the biggest average and the smallest variation
apply(cor(GEBVs[, 3:5]), 2, mean) / apply(cor(GEBVs[, 3:5]), 2, sd)

# defining a threshold to select the individuals to advance
is <- 15 #15%

quantile(GEBVs$SH, ((100 - is) / 100))
GEBVs$Seleted <- GEBVs$SH > quantile(GEBVs$SH, ((100 - is) / 100))
# the number of pre-selected parents
sum(GEBVs$Seleted)
head(GEBVs)
# saving the file
write.csv(GEBVs, "GEBVS.csv")

# 3D plot showing the selected materials and parents
library(ggplot2)
library(ggpubr)
a <- ggplot(GEBVs, aes(x = SDM, 
                       y = SRA,
                       color = Seleted)) + geom_point()
a
#################### the end ####################################
