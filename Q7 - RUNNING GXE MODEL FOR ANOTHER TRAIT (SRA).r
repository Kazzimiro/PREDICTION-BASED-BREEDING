###################################################
        
#           HOMEWORK - HS 541 Plant Breeding Methods - 7TH QUESTION
#  DISCENTE: JOSÉ ARTUR DE OLIVEIRA CASIMIRO
#  DOCENTES: ROBERTO FRITSCHE-NETO E JULIO CESAR DoVALE

#contato: artur.casimiro@alu.ufc.br

#atualização: 24/08/2025

###################################################

#################################### MET analysis #########################
library(breedR)
pheno <- readRDS("pheno")
head(pheno)

# first, we need to create a col for the GxN interaction
pheno$GN <- paste0(pheno$gid, pheno$N)
head(pheno)

# them, get the number of replicates and N levels 
blocks <- length(unique(pheno$rep))
loc <- length(unique(pheno$N))

# Fitting genotype by environment models - with a common variance (diagonal model)
sol <- remlf90(fixed = SRA ~ N + rep, 
                random = ~ gid + GN,
                method = "em",
                data = pheno) 

# heritability and importance of GxE
sol$var
                                        ## FOR SDM ##
#       Estimated variances
##  gid                 0.031630
##  GN                  0.009896
##  Residual            0.020650

                                        ## FOR SRA ##
#       Estimated variances
#gid                5.397e-05
#GN                 1.468e-05
#Residual           1.637e-04



# SDM:
#   A interação Genótipo x Ambiente (GxE) é significativa para a característica SDM.
#   A variância do GN (interação GxE) é alta (0.009896), indicando que o desempenho
#   dos genótipos varia consideravelmente dependendo do ambiente.

# SRA:
#   A interação Genótipo x Ambiente (GxE) é baixa para a característica SRA.
#   A variância do GN (interação GxE) é muito baixa (1.468e-05), sugerindo que
#   o desempenho dos genótipos é mais estável e consistente entre os ambientes.


#########################################################################################

(h2g.plot <- sol$var[1,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
(h2g <- sol$var[1,] / (sol$var[1,] + sol$var[2,]/loc + sol$var[3,]/(loc*blocks)))
(Hgxe <- sol$var[2,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
sol$fit$AIC

                              ## FOR SDM ##

#> (h2g.plot <- sol$var[1,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
#[1] 0.5087172
#> (h2g <- sol$var[1,] / (sol$var[1,] + sol$var[2,]/loc + sol$var[3,]/(loc*blocks)))
#[1] 0.7577772
#> (Hgxe <- sol$var[2,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
#[1] 0.1591611
#> sol$fit$AIC
#[1] -60.36732


                              ## FOR SRA ##
#> (h2g.plot <- sol$var[1,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
#[1] 0.2322789
#> (h2g <- sol$var[1,] / (sol$var[1,] + sol$var[2,]/loc + sol$var[3,]/(loc*blocks)))
#[1] 0.5279014
#> (Hgxe <- sol$var[2,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
#[1] 0.06318055
#> sol$fit$AIC
#[1] -1193.409
# GxE deviations

# SDM:
#   - A herdabilidade do SDM é alta (h2g = 0.758), indicando um forte controle genético.
#   - No entanto, a importância da interação GxE (Hgxe = 0.159) também é significativa.
#   - Isso sugere que o desempenho dos genótipos varia bastante em diferentes ambientes,
#     o que torna a seleção mais complexa.

# SRA:
#   - A herdabilidade do SRA é moderada (h2g = 0.528), menor que a do SDM.
#   - A interação GxE (Hgxe = 0.063) é muito baixa, indicando que o ranqueamento
#     dos genótipos é mais estável e consistente entre os ambientes.
#   - O modelo de SRA teve um ajuste muito melhor (AIC = -1193.4) do que o de SDM,
#     o que reforça a menor complexidade da interação GxE para esta característica.



head(sol$ranef$GN[[1]])
                                 # FOR SDM ##
#              value       s.e.
#10ideal  -0.09963937 0.08332451
#10low     0.03771890 0.08332451
#110ideal -0.04251658 0.08332451
#110low    0.06107710 0.08332451
#112ideal  0.13253822 0.08332451
#112low   -0.05257983 0.08332451

                                 ## FOR SRA ##

#             value       s.e.
#10ideal  -0.00079359 0.00361459
#10low     0.00033951 0.00361457
#110ideal -0.00201923 0.00361459
#110low   -0.00004995 0.00361457
#112ideal -0.00078513 0.00361459
#112low    0.00004392 0.00361457


# SDM:
#   - Variância genética (gid) é significativa (0.031630).
#   - A variância da interação GxE (GN) também é alta (0.009896).
#   - Isso indica que a performance dos genótipos de SDM é fortemente influenciada
#     pela interação com o ambiente.

# SRA:
#   - A variância genética (gid) é muito baixa (5.397e-05).
#   - A variância da interação GxE (GN) também é muito baixa (1.468e-05).
#   - Isso sugere que a característica SRA é mais estável e consistente em
#     diferentes ambientes, com pouca variação genética e GxE.


## Including genomics into the model - only the additive as an example
Ga <- readRDS("Ga")
Za <- model.matrix(~ -1 + gid, data = pheno)
colnames(Za) <- gsub("gid", "", colnames(Za), fixed = T)

sol2 <- remlf90(fixed = SRA ~ N + rep, 
                   random = ~ GN,
                   generic = list(GEBV = list(Za, Ga)),
                   method = "em",
                   data = pheno)

# heritability and importance of GxE
sol2$var
(h2g.plot.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
(h2g.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,]/loc + sol2$var[3,]/(loc*blocks)))
(Hgxe.2 <- sol2$var[1,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
sol2$fit$AIC
                            #####   FOR SDM   ########
#> sol2$var
#Estimated variances
#GN                 0.0413200
#GEBV               0.0002542
#Residual           0.0206500
#> (h2g.plot.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
  #[1] 0.004085227
#> (h2g.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,]/loc + sol2$var[3,]/(loc*blocks)))
  #[1] 0.005453812
#> (Hgxe.2 <- sol2$var[1,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
  #[1] 0.6640503
#> sol2$fit$AIC
  #[1] -35.18977

                                ## FOR SRA ##
#> sol2$var
# Estimated variances
#GN                 4.649e-05
#GEBV               2.114e-05
#Residual           1.639e-04
#> (h2g.plot.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
  #[1] 0.09130566
#> (h2g.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,]/loc + sol2$var[3,]/(loc*blocks)))
  #[1] 0.2156373
#> (Hgxe.2 <- sol2$var[1,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
  #[1] 0.2007947
#> sol2$fit$AIC
  #[1] -1191.724




# SDM:
#   - A importância da interação GxE (Hgxe.2 = 0.664) é extremamente alta.
#   - A herdabilidade (`h2g.2` = 0.005) é praticamente nula.
#   - Isso indica que a maior parte da variação na performance dos genótipos
#     é devido à interação com o ambiente, tornando a seleção para esta
#     característica extremamente difícil.

# SRA:
#   - A herdabilidade (`h2g.2` = 0.216) é moderada, o que é um resultado
#     muito mais favorável para a seleção.
#   - A importância da interação GxE (Hgxe.2 = 0.201) é menor que a herdabilidade.
#   - O modelo para SRA (AIC = -1191.7) se ajusta muito melhor do que o para SDM,
#     confirmando que a GxE é um fator de menor complexidade para esta característica.


# GxE deviations
head(sol2$ranef$GN[[1]])


################## FOR SDM ######################
#             value       s.e.
#10ideal  -0.32240699 0.09559231
#10low    -0.09783596 0.09559231
#110ideal -0.02291405 0.09486829
#110low    0.14645429 0.09486829
#112ideal  0.42054130 0.09481860
#112low    0.11788641 0.09481860


################## FOR SRA ##############

#             value       s.e.
#10ideal  -0.00034835 0.00567436
#10low     0.00235958 0.00567437
#110ideal -0.00659095 0.00556155
#110low   -0.00189219 0.00556156
#112ideal -0.00305403 0.00555693
#112low   -0.00107005 0.00555695



# SDM:
#   - Os desvios (`value`) variam muito, de -0.32 a 0.42.
#   - Isso indica uma interação GxE significativa e imprevisível.
#   - O desempenho de um genótipo muda drasticamente entre os ambientes ideal e de baixo desempenho.

# SRA:
#   - Os desvios (`value`) são muito pequenos e próximos de zero, de -0.006 a 0.002.
#   - Isso sugere que a interação GxE é mínima.
#   - O desempenho dos genótipos de SRA é muito mais estável e consistente em diferentes ambientes.




## Including genomics into the model and into the GxE
Zge <- model.matrix(~ -1 + GN, data = pheno)
(N.levels <- diag(2))
colnames(N.levels) <- rownames(N.levels) <- unique(pheno$N)
N.levels
GxN <- kronecker(N.levels, Ga)
GxN[1:6, 1:6]

sol3 <- remlf90(fixed = SRA ~ N + rep, 
                   random = ~ 1,
                   generic = list(Ga = list(Za, Ga),
                                  GN = list(Zge, precision = GxN)),
                   method = "em",
                   data = pheno)

# heritability and importance of GxE
sol3$var
(h2g.plot.3 <- sol3$var[1,] / (sol3$var[1,] + sol3$var[2,] + sol3$var[3,]))
(h2g.3 <- sol3$var[1,] / (sol3$var[1,] + sol3$var[2,]/loc + sol3$var[3,]/(loc*blocks)))
(Hgxe.3 <- sol3$var[2,] / (sol3$var[1,] + sol3$var[2,] + sol3$var[3,]))
sol3$fit$AIC

# GxE deviations
head(sol3$ranef$GN[[1]])

################################### ERM kernels ################################ 
# If you want to include ERM and GRM into your model
#GxE <- kronecker(ERM, Ga)


####################### Stability and adaptability - Finlay-Wilkinson #########################
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)
library(statgenGxE)

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = pheno, 
                                genotype = "gid", 
                                trial = "N")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, 
                 trait = "SRA")
summary(dropsFW)
 
######################### FOR SDM ###############################
#Environmental effects 
#=====================
#    Trial   EnvEff     SE_EnvEff      EnvMean    SE_EnvMean    Rank
#1   ideal  0.20535    0.007216702   0.8080545    0.1027234      1
#2   low   -0.20535    0.007216702   0.3973545    0.1027234      2

#Anova 
#=====
#                Df Sum  Sq   MeanSq   F value     Pr(>F)    
#Trial           1   9.2771  9.2771   439.5851  < 2.2e-16 ***
#Genotype       54   9.0171  0.1670   7.9123   < 2.2e-16 ***
#Sensitivities  54   2.1841  0.0404   1.9165   0.002066 ** 
#Residual      110   2.3215  0.0211                       
#Total         219   22.7997  0.1041                       
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Most sensitive genotypes
#========================
#      Genotype     GenMean     SE_GenMeaN Rank          Sens     SE_Sens    MSdeviation
#           510     1.02850     0.07263645    1       2.108595  0.3537202  0.02976250
#           312     0.96825     0.07263645    2      1.924763  0.3537202  0.04566425
#             5     0.94000     0.07263645    3      1.921110  0.3537202  0.07054650
#           514     0.84300     0.07263645    4      1.845629  0.3537202  0.03408100
#           511     0.85750     0.07263645    5      1.775018  0.3537202  0.03033000




######################### FOR SRA ###############################



#Environmental effects 
#=====================
#  Trial       EnvEff   SE_EnvEff   EnvMean  SE_EnvMean Rank
#1 ideal  0.003695122 0.000337569 0.1223273 0.009195645    1
#2   low -0.003695122 0.000337569 0.1149364 0.009126908    2

#Anova 
#=====
#                Df Sum Sq    Mean Sq     F value        Pr(>F)    
#Trial           1 0.002986 0.00298566   17.9210      4.822e-05 ***
#Genotype       54 0.022138 0.00040997   2.4608       3.485e-05 ***
#Sensitivities  54 0.010363 0.00019191   1.1519       0.2643    
#Residual      109 0.018159 0.00016660                      
#Total         218 0.053647 0.00024609                      
#--- 
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Most sensitive genotypes
#========================
#  Genotype  GenMean   SE_GenMean   Rank     Sens     SE_Sens   MSdeviation
      #612  0.12625   0.006453699    1   5.344884   1.746546   0.00007025
      #312  0.13800  0.006453699    2   5.277228   1.746546   0.00018450
      #710  0.12875  0.006453699    3   4.668317   1.746546   0.00069725
      #511  0.13475  0.006453699    4   3.585809   1.746546   0.00022525
      #510  0.11800  0.006453699    5   3.382838   1.746546   0.00001250


# SDM:
#   - O efeito ambiental é significativo (Pr < 2.2e-16), com uma grande diferença
#     entre a média dos ambientes 'ideal' e 'low'.
#   - A interação Genótipo x Ambiente (Sensitivities) também é significativa (Pr = 0.002).
#   - Isso confirma que a resposta dos genótipos ao ambiente varia muito,
#     tornando-o um fator de seleção crítico para SDM.

# SRA:
#   - O efeito ambiental é significativo, mas a diferença entre as médias
#     é muito pequena.
#   - A interação GxE (Sensitivities) NÃO é significativa (Pr = 0.2643).
#   - Isso confirma que a resposta dos genótipos ao ambiente é consistente e
#     que a GxE não é um fator importante para a característica SRA.




# let's take a look at the output
names(dropsFW)
dropsFW$estimates
dropsFW$envEffs
dropsFW$fittedGeno



## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")

                    ##################### FINLAY WILKINSON FOR SDM ######################


#Interação Genótipo x Ambiente: Existe uma interação, pois nem todos os genótipos respondem da mesma forma ao ambiente.

#Genótipos Instáveis: A maioria dos genótipos se sai melhor no ambiente ideal, como esperado. As linhas que são mais inclinadas (quase verticais) indicam genótipos mais instáveis.

#Genótipos Estáveis: As linhas que são mais planas mostram genótipos mais estáveis, que têm um desempenho similar nos dois ambientes.

#Genótipos de Alto Desempenho: Há alguns genótipos que se destacam com os valores mais altos de SDM no ambiente ideal.

                   
                   ##################### FINLAY WILKINSON FOR SRA ######################

#Baixa Interação GxE: Ao contrário do SDM, as linhas neste gráfico são muito mais próximas e menos inclinadas. Isso sugere que a interação genótipo x ambiente é baixa para a SRA. O desempenho dos genótipos é relativamente consistente entre os dois ambientes.

#Genótipos Mais Estáveis: A maioria dos genótipos demonstra uma estabilidade maior. As linhas que são mais horizontais (e menos inclinadas) indicam genótipos que mantêm um desempenho similar, independentemente do ambiente.

#Seleção: A seleção para a SRA seria mais simples, pois o desempenho de um genótipo em um ambiente de baixa performance tende a ser um bom indicador de seu desempenho em um ambiente ideal.





############################# GGE-Biplot Analysis ######################################
library(metan)

model.gge <- gge(pheno, N, gid, SRA, svp = "symmetrical")

model.gge$SDM

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where

######## the end ##########
