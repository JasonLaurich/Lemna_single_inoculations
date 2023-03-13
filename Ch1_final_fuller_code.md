---
title: "Final code with more full exploration of data and errors"
author: "JasonL"
date: "2023-02-06"
output: 
  html_document: 
    keep_md: yes
---

SECTION 1: Upload packages


```r
library(lme4) # mixed effects models
library(lmerTest) # get results of lmer models
library(reshape2) # melting dataframes
library(plyr) # manipulating data frames
library(dplyr) # manipulating data frames
library(ggplot2) # plotting
library(cowplot) # luxury plotting
library(grid) # arranging and labelling plots
library(gridExtra)
library(emmeans) # allows extraction of bacteria-specific regression estimates.
library(ggrepel)
```

SECTION 2: Upload data for 2020 and 2019 experiments.


```r
setwd("C:/Users/jason/Dropbox/My PC (DESKTOP-IOO274E)/Desktop/Data/Lemna_singleinoculations/Lemna_single_inoculations")

df.CW<-read.table('CW_mixed_inoculation_data.txt', header=F, na.strings="", sep="\t")
names(df.CW)<-c('pop','plt','bac','id','frd0','pix0','per0','grn0','grnper0','frd1','pix1','per1','grn1','grnper1','abs')

#Fix the structure for bacterial inoculation and plant
df.CW$bac <- factor(df.CW$bac, levels = c("Control", "1", "2", "3", "4", "5","6","7","8","9","10","All10"))
df.CW$plt <- as.factor(df.CW$plt)

# Add in edge effects. So for each 24 well plate, all but 6,7,10,11,14,15,18, and 19 are edge.
edge<-c(rep("Y",5), rep("N",2), rep("Y",2),rep("N",2), rep("Y",2), rep("N",2), rep("Y",2),rep("N",2), rep("Y",5))

#Final plate has 8 samples only
edge2<-rep(edge,38)
edge2<-c(edge2, rep("Y",5), rep("N",2), "Y")
df.CW$edge<-edge2

df.CW$edge <- as.factor(df.CW$edge)

plate<-vector()
# Create vector for plate - 38 full, 39th with 8 wells
for (i in 1:38){
  plate<-c(plate, rep(i,24))
}
plate<-c(plate,rep(39,8))

df.CW$plate <- plate
df.CW$plate <- as.factor(df.CW$plate)

#Take out errors (pipetting errors, wrong bacteria)
df.CW<-df.CW[-c(453,494,798),]

# Calculate growth.
df.CW$Grt<-df.CW$pix1 - df.CW$pix0

# Let's log abs count
df.CW$logabs<-log(df.CW$abs)

#Let's create separate dataframes for each population # Going to get the label 2020
df.C.2020<-subset(df.CW, df.CW$pop== "Churchill")
df.W.2020<-subset(df.CW, df.CW$pop== "Wellspring")

# Change factor order of bacteria to reflect their effects on duckweed growth.
#First, let's change the order of bacteria, ranked by effect. 
df.C.2020$bac<- factor(df.C.2020$bac, levels = c("Control","1","8","6","3","7","10","9","4","5","2","All10"))
# Assign names from sequencing
levels(df.C.2020$bac) <- c("Control","Flavobacterium \nsuccinicans 1","Bosea \nmassiliensis","Aeromonas \nsalmonicida","Ohtaekwangia \nkoreensis","Flavobacterium \nsuccinicans 2","Falsiroseomonas \nstagni","Parasediminibacterium \npaludis","Arcicella sp.","Microbacterium \noxydans","Pseudomonas \nprotogens","All 10 bacteria")

df.W.2020$bac<- factor(df.W.2020$bac, levels = c("Control","2","8","1","9","6","3","5","4","7","10","All10"))
levels(df.W.2020$bac) <- c("Control","Sphingomonas \npituitosa 1","Flaviflagellibacter \ndeserti","Rhizobium \nrosettiformans","Rhizorhabdus \nwittichii 1","Rhizobium \ncapsici 1","Pseudomonas \nprotogens","Rhizorhabdus \nwittichii 2","Rhizobium \ncapsici 2","Sphingomonas \npituitosa 2","Fervidobacterium \nriparium","All 10 bacteria")

#####################Load 2019 data#######################

df.CKW<-read.table('July2019CKW.txt', header=F, na.strings="", sep="\t")
names(df.CKW)<-c('pop','bac','id','frd0','pix0','per0','grn0','grnper0','frd1','pix1','per1','grn1','grnper1','notes','algae_lvl')

#Fix the structure
df.CKW$pop<-as.factor(df.CKW$pop)
df.CKW$bac<-as.factor(df.CKW$bac)

# Add in edge (15 plates)
edge3<-rep(edge,15)
df.CKW$edge<-edge3

df.CKW$edge <- as.factor(df.CKW$edge)
df.CKW$algae_lvl<-as.factor(df.CKW$algae_lvl)

plate2<-vector()
# Create vector for plate
for (i in 1:15){
  plate2<-c(plate2, rep(i,24))
}

df.CKW$plate <- plate2
df.CKW$plate <- as.factor(df.CKW$plate)

df.CKW$Grt<-df.CKW$pix1 - df.CKW$pix0

# We want to include algal level as an ordinal variable.
df.CKW$alg.ord<- factor(df.CKW$algae_lvl, order = TRUE, levels = c("None","Low", "Med", "High"))

# Let's break this down by population.
df.C.2019<-subset(df.CKW, df.CKW$pop== "Churchill")
df.K.2019<-subset(df.CKW, df.CKW$pop== "Kelso")
df.W.2019<-subset(df.CKW, df.CKW$pop== "Wellspring ")

# Change factor order of bacteria to reflect their effects on duckweed growth.
#First, let's change the order of bacteria, ranked by effect. 
df.C.2019$bac<- factor(df.C.2019$bac, levels = c("0","7","8","5","6","10","4","9","2","1","3","All"))
# Assign names from sequencing
levels(df.C.2019$bac) <- c("Control","Aeromonas \nsalmonicida 1","Aeromonas \nsalmonicida 2","Microbacterium sp.","Parasediminibacterium \npaludis","Flavibacterium \nsuccinicans 1","Arcicella sp.","Falsiroseomonas \nstagni","Pseudomonas \nprotogens","Flavibacterium \nsuccinicans 2","Ohtaekwangia \nkoreensis","All 10 bacteria")

df.K.2019$bac<- factor(df.K.2019$bac, levels = c("0","8","6","7","3","9","2","5","1","4","10","All"))
levels(df.K.2019$bac) <- c("Control","Pseudomonas \nprotogens","Flavobacterium \nbranchiicola","Pseudomonas \nbaetica","Rhizobium sp.","Pseudacidovorax \nintermedius","Rhizobium \nrosettiformans","Flavobacterium \npanici 1","Flavobacterium \npanici 2","Flavobacterium \npanici 3","Pseudomonas \nsesami","All 10 bacteria")

df.W.2019$bac<- factor(df.W.2019$bac, levels = c("0","8","5","3","2","4","6","9","7","1","10","All"))
levels(df.W.2019$bac) <- c("Control","Rhizobium sp.","Agrobacterium sp.","Rhizobium \nrosettiformans 1","Rhizobium \nrhizophilum","Rhizobium \nrosettiformans 2","Rhizobium \nhelianthi","Sphingomonas \npituitosa","Rhizorhabdus \nwittichii","Rhizobium \nglycinendophyticum","Phenylobacterium \npanacis","All 10 bacteria")
```

SECTION 3: Plant growth models and figures


```r
#Churchill 2019

lmer.grt.C2019<-lmer(pix1 ~ pix0 + bac + alg.ord + (1|edge) + (1|plate), data=df.C.2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.C2019)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##             Sum Sq    Mean Sq NumDF   DenDF F value   Pr(>F)    
## pix0    1080020412 1080020412     1 103.388 80.8212 1.26e-14 ***
## bac      100512555    9137505    11  93.169  0.6838   0.7506    
## alg.ord      58434      58434     1 104.601  0.0044   0.9474    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.C2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate)
##             npar  logLik    AIC    LRT Df Pr(>Chisq)    
## <none>        17 -1047.0 2127.9                         
## (1 | edge)    16 -1052.1 2136.2 10.288  1  0.0013388 ** 
## (1 | plate)   16 -1053.2 2138.3 12.378  1  0.0004343 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.C2019)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate)
##    Data: df.C.2019
## 
## REML criterion at convergence: 2093.9
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.5423 -0.5087 -0.1239  0.4692  2.8763 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept)  6439616 2538    
##  edge     (Intercept)  4286331 2070    
##  Residual             13363079 3656    
## Number of obs: 119, groups:  plate, 15; edge, 2
## 
## Fixed effects:
##                                      Estimate Std. Error         df t value
## (Intercept)                         3476.3537  2200.3682     4.3688   1.580
## pix0                                   2.3243     0.2585   103.3878   8.990
## bacAeromonas \nsalmonicida 1       -1171.2306  1772.6294    96.1508  -0.661
## bacAeromonas \nsalmonicida 2        1045.0554  1698.0783    92.2427   0.615
## bacMicrobacterium sp.               -253.5136  1793.0124    97.4365  -0.141
## bacParasediminibacterium \npaludis   -83.2679  1696.7590    92.1917  -0.049
## bacFlavibacterium \nsuccinicans 1    396.5644  1714.7060    93.4988   0.231
## bacArcicella sp.                    -410.7178  1797.0384    97.8940  -0.229
## bacFalsiroseomonas \nstagni          384.4008  1752.6485    92.6260   0.219
## bacPseudomonas \nprotogens          1168.3567  1710.2152    93.1562   0.683
## bacFlavibacterium \nsuccinicans 2   1561.2487  1710.6465    93.0179   0.913
## bacOhtaekwangia \nkoreensis         2570.7462  1738.1356    94.1128   1.479
## bacAll 10 bacteria                   557.9978  1731.8132    93.5047   0.322
## alg.ord.L                            -54.1890   819.4658   104.6014  -0.066
##                                    Pr(>|t|)    
## (Intercept)                           0.183    
## pix0                               1.26e-14 ***
## bacAeromonas \nsalmonicida 1          0.510    
## bacAeromonas \nsalmonicida 2          0.540    
## bacMicrobacterium sp.                 0.888    
## bacParasediminibacterium \npaludis    0.961    
## bacFlavibacterium \nsuccinicans 1     0.818    
## bacArcicella sp.                      0.820    
## bacFalsiroseomonas \nstagni           0.827    
## bacPseudomonas \nprotogens            0.496    
## bacFlavibacterium \nsuccinicans 2     0.364    
## bacOhtaekwangia \nkoreensis           0.142    
## bacAll 10 bacteria                    0.748    
## alg.ord.L                             0.947    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 14 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.C2019)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-1.png)<!-- -->

```r
# Try running an ANCOVA style model (e.g. checking for variaiton in slope effects of pix0vpix1 for each bac.)
lmer.grt.C2019.2<-lmer(pix1 ~ pix0*bac + alg.ord + (1|edge) + (1|plate), data=df.C.2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.C2019.2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##             Sum Sq   Mean Sq NumDF  DenDF F value    Pr(>F)    
## pix0     822787827 822787827     1 91.439 60.3437 1.125e-11 ***
## bac       72248785   6568071    11 82.876  0.4817    0.9097    
## alg.ord     607852    607852     1 93.489  0.0446    0.8332    
## pix0:bac 113356423  10305129    11 83.607  0.7558    0.6822    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.C2019.2)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate) + pix0:bac
##             npar  logLik    AIC     LRT Df Pr(>Chisq)    
## <none>        28 -1033.8 2123.7                          
## (1 | edge)    27 -1038.7 2131.3  9.6503  1  0.0018932 ** 
## (1 | plate)   27 -1039.7 2133.5 11.8163  1  0.0005871 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.C2019.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 * bac + alg.ord + (1 | edge) + (1 | plate)
##    Data: df.C.2019
## 
## REML criterion at convergence: 2067.7
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.12448 -0.56640 -0.07159  0.46383  2.59389 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept)  6924341 2631    
##  edge     (Intercept)  4674342 2162    
##  Residual             13635029 3693    
## Number of obs: 119, groups:  plate, 15; edge, 2
## 
## Fixed effects:
##                                          Estimate Std. Error        df t value
## (Intercept)                             1594.5781  3153.6527   14.2384   0.506
## pix0                                       2.9415     0.7519   87.8231   3.912
## bacAeromonas \nsalmonicida 1            2701.6705  3897.8636   82.1784   0.693
## bacAeromonas \nsalmonicida 2            2226.5438  3359.6507   83.7103   0.663
## bacMicrobacterium sp.                   1411.2225  3952.1772   85.1540   0.357
## bacParasediminibacterium \npaludis      5572.0499  3770.9113   81.2634   1.478
## bacFlavibacterium \nsuccinicans 1       5860.9443  4157.0619   80.5403   1.410
## bacArcicella sp.                         670.9082  3882.8848   87.7352   0.173
## bacFalsiroseomonas \nstagni             2136.9536  3666.2409   83.4749   0.583
## bacPseudomonas \nprotogens              1738.6842  3506.9549   82.7515   0.496
## bacFlavibacterium \nsuccinicans 2        213.1152  3579.5040   81.4865   0.060
## bacOhtaekwangia \nkoreensis             3804.1816  3384.3778   81.1921   1.124
## bacAll 10 bacteria                      3106.1972  3415.0501   81.3158   0.910
## alg.ord.L                               -183.2724   868.0130   93.4888  -0.211
## pix0:bacAeromonas \nsalmonicida 1         -1.3605     1.2560   82.3908  -1.083
## pix0:bacAeromonas \nsalmonicida 2         -0.3464     0.9305   81.4381  -0.372
## pix0:bacMicrobacterium sp.                -0.5649     1.4313   84.8937  -0.395
## pix0:bacParasediminibacterium \npaludis   -1.7607     1.0505   84.0939  -1.676
## pix0:bacFlavibacterium \nsuccinicans 1    -1.9073     1.2962   80.9892  -1.471
## pix0:bacArcicella sp.                     -0.3942     1.0156   88.6489  -0.388
## pix0:bacFalsiroseomonas \nstagni          -0.6168     1.0409   82.0390  -0.593
## pix0:bacPseudomonas \nprotogens           -0.1889     1.0061   83.7515  -0.188
## pix0:bacFlavibacterium \nsuccinicans 2     0.7210     1.1756   83.3460   0.613
## pix0:bacOhtaekwangia \nkoreensis          -0.3528     1.0951   81.2377  -0.322
## pix0:bacAll 10 bacteria                   -0.8811     0.9289   81.1756  -0.949
##                                         Pr(>|t|)    
## (Intercept)                              0.62085    
## pix0                                     0.00018 ***
## bacAeromonas \nsalmonicida 1             0.49019    
## bacAeromonas \nsalmonicida 2             0.50932    
## bacMicrobacterium sp.                    0.72192    
## bacParasediminibacterium \npaludis       0.14337    
## bacFlavibacterium \nsuccinicans 1        0.16243    
## bacArcicella sp.                         0.86322    
## bacFalsiroseomonas \nstagni              0.56155    
## bacPseudomonas \nprotogens               0.62136    
## bacFlavibacterium \nsuccinicans 2        0.95267    
## bacOhtaekwangia \nkoreensis              0.26431    
## bacAll 10 bacteria                       0.36574    
## alg.ord.L                                0.83324    
## pix0:bacAeromonas \nsalmonicida 1        0.28187    
## pix0:bacAeromonas \nsalmonicida 2        0.71067    
## pix0:bacMicrobacterium sp.               0.69406    
## pix0:bacParasediminibacterium \npaludis  0.09743 .  
## pix0:bacFlavibacterium \nsuccinicans 1   0.14504    
## pix0:bacArcicella sp.                    0.69886    
## pix0:bacFalsiroseomonas \nstagni         0.55513    
## pix0:bacPseudomonas \nprotogens          0.85150    
## pix0:bacFlavibacterium \nsuccinicans 2   0.54136    
## pix0:bacOhtaekwangia \nkoreensis         0.74814    
## pix0:bacAll 10 bacteria                  0.34566    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 25 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.C2019.2)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-2.png)<!-- -->

```r
#Kelso 2019

lmer.grt.K2019<-lmer(pix1 ~ pix0 + bac + alg.ord + (1|edge) + (1|plate), data=df.K.2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.K2019)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##             Sum Sq    Mean Sq NumDF  DenDF  F value    Pr(>F)    
## pix0    1.0097e+10 1.0097e+10     1 94.637 142.7559 < 2.2e-16 ***
## bac     3.0854e+09 2.8049e+08    11 96.376   3.9659 9.153e-05 ***
## alg.ord 5.7150e+08 2.8575e+08     2 31.719   4.0402   0.02735 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.K2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate)
##             npar  logLik    AIC     LRT Df Pr(>Chisq)    
## <none>        18 -1123.0 2282.0                          
## (1 | edge)    17 -1129.5 2293.0 13.0512  1  0.0003031 ***
## (1 | plate)   17 -1127.3 2288.7  8.7322  1  0.0031264 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.K2019)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate)
##    Data: df.K.2019
## 
## REML criterion at convergence: 2246
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.9699 -0.5435 -0.1227  0.5151  3.1955 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept) 19968963 4469    
##  edge     (Intercept) 32543548 5705    
##  Residual             70726146 8410    
## Number of obs: 119, groups:  plate, 15; edge, 2
## 
## Fixed effects:
##                                    Estimate Std. Error         df t value
## (Intercept)                      -2564.5047  5262.5079     2.5668  -0.487
## pix0                                 3.5913     0.3006    94.6371  11.948
## bacPseudomonas \nprotogens       -9601.0933  3852.7821    93.8857  -2.492
## bacFlavobacterium \nbranchiicola -5598.8372  4149.9362    96.3171  -1.349
## bacPseudomonas \nbaetica          -667.3961  4088.1359    99.6765  -0.163
## bacRhizobium sp.                  -972.8509  3898.6045    94.7027  -0.250
## bacPseudacidovorax \nintermedius   111.7976  3946.3040    95.6852   0.028
## bacRhizobium \nrosettiformans     4340.8435  3899.6644    96.6285   1.113
## bacFlavobacterium \npanici 1      2348.2983  3985.6073    97.2631   0.589
## bacFlavobacterium \npanici 2      8036.7149  3956.5186    96.1191   2.031
## bacFlavobacterium \npanici 3      6689.4345  4042.2252    99.1663   1.655
## bacPseudomonas \nsesami           7443.1236  4005.0536    95.4206   1.858
## bacAll 10 bacteria               -6051.1293  3954.0625    98.9111  -1.530
## alg.ord.L                        -7617.2145  2693.3356    30.3344  -2.828
## alg.ord.Q                          409.4683  2235.5886    95.4505   0.183
##                                  Pr(>|t|)    
## (Intercept)                       0.66456    
## pix0                              < 2e-16 ***
## bacPseudomonas \nprotogens        0.01446 *  
## bacFlavobacterium \nbranchiicola  0.18046    
## bacPseudomonas \nbaetica          0.87065    
## bacRhizobium sp.                  0.80348    
## bacPseudacidovorax \nintermedius  0.97746    
## bacRhizobium \nrosettiformans     0.26841    
## bacFlavobacterium \npanici 1      0.55710    
## bacFlavobacterium \npanici 2      0.04499 *  
## bacFlavobacterium \npanici 3      0.10111    
## bacPseudomonas \nsesami           0.06619 .  
## bacAll 10 bacteria                0.12912    
## alg.ord.L                         0.00822 ** 
## alg.ord.Q                         0.85506    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 15 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.K2019)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-3.png)<!-- -->

```r
lmer.grt.K2019.2<-lmer(pix1 ~ pix0*bac + alg.ord + (1|edge) + (1|plate), data=df.K.2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.K2019.2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##              Sum Sq    Mean Sq NumDF  DenDF  F value    Pr(>F)    
## pix0     7196250397 7196250397     1 77.861 132.4832 < 2.2e-16 ***
## bac       837939168   76176288    11 88.611   1.4024  0.185918    
## alg.ord   626961301  313480650     2 27.748   5.7712  0.008022 ** 
## pix0:bac 2442992316  222090211    11 88.160   4.0887 7.495e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.K2019.2)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate) + pix0:bac
##             npar  logLik    AIC    LRT Df Pr(>Chisq)   
## <none>        29 -1093.2 2244.3                        
## (1 | edge)    28 -1097.5 2250.9 8.6130  1   0.003338 **
## (1 | plate)   28 -1095.8 2247.7 5.3942  1   0.020203 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.K2019.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 * bac + alg.ord + (1 | edge) + (1 | plate)
##    Data: df.K.2019
## 
## REML criterion at convergence: 2186.3
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.9238 -0.5074 -0.0984  0.3675  3.7799 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept) 12647588 3556    
##  edge     (Intercept) 18820342 4338    
##  Residual             54318198 7370    
## Number of obs: 119, groups:  plate, 15; edge, 2
## 
## Fixed effects:
##                                         Estimate Std. Error         df t value
## (Intercept)                            3.975e+03  5.341e+03  7.525e+00   0.744
## pix0                                   2.176e+00  8.148e-01  9.085e+01   2.670
## bacPseudomonas \nprotogens            -7.251e+03  5.435e+03  8.439e+01  -1.334
## bacFlavobacterium \nbranchiicola      -3.103e+03  5.790e+03  9.063e+01  -0.536
## bacPseudomonas \nbaetica              -1.259e+04  5.552e+03  9.192e+01  -2.268
## bacRhizobium sp.                      -5.365e+03  5.597e+03  8.718e+01  -0.959
## bacPseudacidovorax \nintermedius      -1.339e+04  5.945e+03  8.935e+01  -2.253
## bacRhizobium \nrosettiformans          1.116e+03  5.728e+03  8.744e+01   0.195
## bacFlavobacterium \npanici 1          -6.079e+03  5.483e+03  8.764e+01  -1.109
## bacFlavobacterium \npanici 2           2.829e+03  5.994e+03  8.524e+01   0.472
## bacFlavobacterium \npanici 3          -5.532e+03  7.285e+03  8.602e+01  -0.759
## bacPseudomonas \nsesami               -6.878e+03  5.487e+03  8.918e+01  -1.253
## bacAll 10 bacteria                    -7.425e+03  5.550e+03  9.181e+01  -1.338
## alg.ord.L                             -7.920e+03  2.341e+03  2.634e+01  -3.382
## alg.ord.Q                             -7.447e+00  2.006e+03  8.808e+01  -0.004
## pix0:bacPseudomonas \nprotogens       -9.499e-01  1.104e+00  8.169e+01  -0.860
## pix0:bacFlavobacterium \nbranchiicola -2.022e+00  1.418e+00  9.198e+01  -1.426
## pix0:bacPseudomonas \nbaetica          3.567e+00  1.389e+00  9.090e+01   2.568
## pix0:bacRhizobium sp.                  1.048e+00  1.059e+00  8.738e+01   0.989
## pix0:bacPseudacidovorax \nintermedius  3.451e+00  1.250e+00  8.921e+01   2.761
## pix0:bacRhizobium \nrosettiformans     8.417e-01  1.192e+00  8.561e+01   0.706
## pix0:bacFlavobacterium \npanici 1      1.937e+00  9.955e-01  8.591e+01   1.946
## pix0:bacFlavobacterium \npanici 2      1.177e+00  1.177e+00  8.621e+01   1.000
## pix0:bacFlavobacterium \npanici 3      3.680e+00  1.891e+00  8.643e+01   1.946
## pix0:bacPseudomonas \nsesami           3.552e+00  1.037e+00  8.978e+01   3.426
## pix0:bacAll 10 bacteria                4.366e-01  1.013e+00  9.196e+01   0.431
##                                       Pr(>|t|)    
## (Intercept)                           0.479288    
## pix0                                  0.008986 ** 
## bacPseudomonas \nprotogens            0.185764    
## bacFlavobacterium \nbranchiicola      0.593390    
## bacPseudomonas \nbaetica              0.025663 *  
## bacRhizobium sp.                      0.340412    
## bacPseudacidovorax \nintermedius      0.026710 *  
## bacRhizobium \nrosettiformans         0.845986    
## bacFlavobacterium \npanici 1          0.270633    
## bacFlavobacterium \npanici 2          0.638211    
## bacFlavobacterium \npanici 3          0.449687    
## bacPseudomonas \nsesami               0.213303    
## bacAll 10 bacteria                    0.184218    
## alg.ord.L                             0.002257 ** 
## alg.ord.Q                             0.997047    
## pix0:bacPseudomonas \nprotogens       0.392243    
## pix0:bacFlavobacterium \nbranchiicola 0.157358    
## pix0:bacPseudomonas \nbaetica         0.011876 *  
## pix0:bacRhizobium sp.                 0.325235    
## pix0:bacPseudacidovorax \nintermedius 0.007003 ** 
## pix0:bacRhizobium \nrosettiformans    0.481933    
## pix0:bacFlavobacterium \npanici 1     0.054915 .  
## pix0:bacFlavobacterium \npanici 2     0.320080    
## pix0:bacFlavobacterium \npanici 3     0.054879 .  
## pix0:bacPseudomonas \nsesami          0.000925 ***
## pix0:bacAll 10 bacteria               0.667574    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 26 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.K2019.2)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-4.png)<!-- -->

```r
#Wellspring 2019

lmer.grt.W2019<-lmer(pix1 ~ pix0 + bac + alg.ord + (1|edge) + (1|plate), data=df.W.2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.W2019)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##            Sum Sq   Mean Sq NumDF  DenDF  F value Pr(>F)    
## pix0    616242532 616242532     1 84.200 141.8587 <2e-16 ***
## bac      75414969   6855906    11 95.486   1.5782 0.1176    
## alg.ord  12621366   4207122     3 84.247   0.9685 0.4116    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.W2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate)
##             npar  logLik    AIC     LRT Df Pr(>Chisq)
## <none>        19 -935.11 1908.2                      
## (1 | edge)    18 -935.43 1906.9 0.64679  1     0.4213
## (1 | plate)   18 -935.51 1907.0 0.80646  1     0.3692
```

```r
summary(lmer.grt.W2019)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate)
##    Data: df.W.2019
## 
## REML criterion at convergence: 1870.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.0295 -0.4799 -0.0527  0.4288  4.1734 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept)  290497   539    
##  edge     (Intercept)  187497   433    
##  Residual             4344058  2084    
## Number of obs: 116, groups:  plate, 15; edge, 2
## 
## Fixed effects:
##                                     Estimate Std. Error         df t value
## (Intercept)                        2535.2167  1083.7966    40.6148   2.339
## pix0                                  1.8151     0.1524    84.2005  11.910
## bacRhizobium sp.                  -1439.2679  1023.6353    99.9616  -1.406
## bacAgrobacterium sp.              -1243.0270  1018.0656    95.0318  -1.221
## bacRhizobium \nrosettiformans 1   -1484.7596  1021.8051    95.4852  -1.453
## bacRhizobium \nrhizophilum         -742.6109   988.1205    96.3971  -0.752
## bacRhizobium \nrosettiformans 2    -747.2002   998.8351    99.6628  -0.748
## bacRhizobium \nhelianthi           -611.5977  1011.5788    96.9725  -0.605
## bacSphingomonas \npituitosa          19.7042   982.0143    96.9734   0.020
## bacRhizorhabdus \nwittichii        -158.6976  1003.8651    99.3694  -0.158
## bacRhizobium \nglycinendophyticum  1101.4995  1014.9069    99.3208   1.085
## bacPhenylobacterium \npanacis      1208.1404   979.6191    95.2702   1.233
## bacAll 10 bacteria                 -199.9646  1045.6015    95.4664  -0.191
## alg.ord.L                          1112.7184  1686.5472    96.1303   0.660
## alg.ord.Q                           785.6508  1325.3524    99.7440   0.593
## alg.ord.C                           -30.4586   647.9390    96.4855  -0.047
##                                   Pr(>|t|)    
## (Intercept)                         0.0243 *  
## pix0                                <2e-16 ***
## bacRhizobium sp.                    0.1628    
## bacAgrobacterium sp.                0.2251    
## bacRhizobium \nrosettiformans 1     0.1495    
## bacRhizobium \nrhizophilum          0.4542    
## bacRhizobium \nrosettiformans 2     0.4562    
## bacRhizobium \nhelianthi            0.5469    
## bacSphingomonas \npituitosa         0.9840    
## bacRhizorhabdus \nwittichii         0.8747    
## bacRhizobium \nglycinendophyticum   0.2804    
## bacPhenylobacterium \npanacis       0.2205    
## bacAll 10 bacteria                  0.8487    
## alg.ord.L                           0.5110    
## alg.ord.Q                           0.5547    
## alg.ord.C                           0.9626    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 16 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.W2019)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-5.png)<!-- -->

```r
lmer.grt.W2019.2<-lmer(pix1 ~ pix0*bac + alg.ord + (1|edge) + (1|plate), data=df.W.2019)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.W2019.2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##             Sum Sq   Mean Sq NumDF  DenDF F value    Pr(>F)    
## pix0     307926948 307926948     1 80.831 75.3647 3.444e-13 ***
## bac       50292184   4572017    11 85.712  1.1190   0.35681    
## alg.ord   30890952  10296984     3 48.394  2.5202   0.06890 .  
## pix0:bac  79899281   7263571    11 85.697  1.7777   0.07043 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.W2019.2)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + alg.ord + (1 | edge) + (1 | plate) + pix0:bac
##             npar logLik    AIC     LRT Df Pr(>Chisq)
## <none>        30 -920.8 1901.6                      
## (1 | edge)    29 -920.8 1899.6 0.00011  1     0.9916
## (1 | plate)   29 -921.1 1900.2 0.59070  1     0.4421
```

```r
summary(lmer.grt.W2019.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 * bac + alg.ord + (1 | edge) + (1 | plate)
##    Data: df.W.2019
## 
## REML criterion at convergence: 1841.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.1285 -0.5095 -0.0343  0.3538  3.9652 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept)  240083   489.98 
##  edge     (Intercept)    1875    43.31 
##  Residual             4085825  2021.34 
## Number of obs: 116, groups:  plate, 15; edge, 2
## 
## Fixed effects:
##                                          Estimate Std. Error         df t value
## (Intercept)                             3.443e+03  2.049e+03  8.847e+01   1.680
## pix0                                    1.531e+00  6.602e-01  8.567e+01   2.319
## bacRhizobium sp.                        1.104e+03  2.824e+03  8.654e+01   0.391
## bacAgrobacterium sp.                   -2.348e+03  2.482e+03  8.668e+01  -0.946
## bacRhizobium \nrosettiformans 1         1.416e+03  2.424e+03  8.888e+01   0.584
## bacRhizobium \nrhizophilum             -1.942e+03  2.170e+03  8.764e+01  -0.895
## bacRhizobium \nrosettiformans 2        -1.586e+03  2.365e+03  8.840e+01  -0.671
## bacRhizobium \nhelianthi               -2.150e+03  2.253e+03  8.584e+01  -0.954
## bacSphingomonas \npituitosa            -1.775e+03  2.301e+03  8.434e+01  -0.772
## bacRhizorhabdus \nwittichii            -1.486e+03  2.285e+03  8.779e+01  -0.650
## bacRhizobium \nglycinendophyticum      -3.344e+03  2.548e+03  8.810e+01  -1.313
## bacPhenylobacterium \npanacis           1.678e+03  2.538e+03  8.837e+01   0.661
## bacAll 10 bacteria                     -7.118e+02  2.217e+03  8.718e+01  -0.321
## alg.ord.L                               1.163e+03  1.633e+03  8.545e+01   0.712
## alg.ord.Q                               1.020e+03  1.273e+03  8.340e+01   0.802
## alg.ord.C                              -3.726e+02  6.303e+02  8.049e+01  -0.591
## pix0:bacRhizobium sp.                  -1.914e+00  1.511e+00  8.623e+01  -1.267
## pix0:bacAgrobacterium sp.               5.279e-01  9.859e-01  8.541e+01   0.535
## pix0:bacRhizobium \nrosettiformans 1   -1.218e+00  8.573e-01  8.803e+01  -1.421
## pix0:bacRhizobium \nrhizophilum         4.934e-01  7.805e-01  8.712e+01   0.632
## pix0:bacRhizobium \nrosettiformans 2    3.339e-01  8.271e-01  8.613e+01   0.404
## pix0:bacRhizobium \nhelianthi           5.631e-01  7.346e-01  8.272e+01   0.766
## pix0:bacSphingomonas \npituitosa        6.516e-01  7.717e-01  8.393e+01   0.844
## pix0:bacRhizorhabdus \nwittichii        4.865e-01  7.431e-01  8.757e+01   0.655
## pix0:bacRhizobium \nglycinendophyticum  2.808e+00  1.236e+00  8.717e+01   2.272
## pix0:bacPhenylobacterium \npanacis     -2.592e-01  8.846e-01  8.709e+01  -0.293
## pix0:bacAll 10 bacteria                 9.629e-02  7.855e-01  8.693e+01   0.123
##                                        Pr(>|t|)  
## (Intercept)                              0.0964 .
## pix0                                     0.0228 *
## bacRhizobium sp.                         0.6967  
## bacAgrobacterium sp.                     0.3468  
## bacRhizobium \nrosettiformans 1          0.5606  
## bacRhizobium \nrhizophilum               0.3733  
## bacRhizobium \nrosettiformans 2          0.5042  
## bacRhizobium \nhelianthi                 0.3427  
## bacSphingomonas \npituitosa              0.4425  
## bacRhizorhabdus \nwittichii              0.5172  
## bacRhizobium \nglycinendophyticum        0.1927  
## bacPhenylobacterium \npanacis            0.5103  
## bacAll 10 bacteria                       0.7489  
## alg.ord.L                                0.4781  
## alg.ord.Q                                0.4250  
## alg.ord.C                                0.5561  
## pix0:bacRhizobium sp.                    0.2085  
## pix0:bacAgrobacterium sp.                0.5937  
## pix0:bacRhizobium \nrosettiformans 1     0.1589  
## pix0:bacRhizobium \nrhizophilum          0.5290  
## pix0:bacRhizobium \nrosettiformans 2     0.6875  
## pix0:bacRhizobium \nhelianthi            0.4456  
## pix0:bacSphingomonas \npituitosa         0.4009  
## pix0:bacRhizorhabdus \nwittichii         0.5144  
## pix0:bacRhizobium \nglycinendophyticum   0.0255 *
## pix0:bacPhenylobacterium \npanacis       0.7702  
## pix0:bacAll 10 bacteria                  0.9027  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 27 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.W2019.2)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-6.png)<!-- -->

```r
#Churchill 2020

lmer.grt.C2020<-lmer(pix1 ~ pix0 + bac + (1|edge) + (1|plate), data=df.C.2020)
anova(lmer.grt.C2020)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##          Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)    
## pix0 4295384852 4295384852     1 220.01 86.4526 < 2.2e-16 ***
## bac  1353542373  123049307    11 207.76  2.4766  0.006176 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.C2020)
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + (1 | edge) + (1 | plate)
##             npar  logLik    AIC     LRT Df Pr(>Chisq)    
## <none>        16 -2345.2 4722.4                          
## (1 | edge)    15 -2350.9 4731.8 11.4585  1  0.0007117 ***
## (1 | plate)   15 -2349.6 4729.1  8.7383  1  0.0031159 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.C2020)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 + bac + (1 | edge) + (1 | plate)
##    Data: df.C.2020
## 
## REML criterion at convergence: 4690.4
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.56199 -0.61502  0.02175  0.69973  3.09833 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept)  9462790 3076    
##  edge     (Intercept)  7800174 2793    
##  Residual             49684859 7049    
## Number of obs: 237, groups:  plate, 39; edge, 2
## 
## Fixed effects:
##                                     Estimate Std. Error        df t value
## (Intercept)                        3878.8193  3030.8194    4.8295   1.280
## pix0                                  4.9076     0.5278  220.0131   9.298
## bacFlavobacterium \nsuccinicans 1  -313.0947  2380.2982  211.3749  -0.132
## bacBosea \nmassiliensis            1614.9724  2355.4880  216.7209   0.686
## bacAeromonas \nsalmonicida         -182.9379  2352.5254  217.5792  -0.078
## bacOhtaekwangia \nkoreensis        -558.6037  2337.8576  213.4642  -0.239
## bacFlavobacterium \nsuccinicans 2  1814.4222  2357.8891  216.2553   0.770
## bacFalsiroseomonas \nstagni         707.6477  2310.0205  207.2307   0.306
## bacParasediminibacterium \npaludis 2265.6999  2314.7662  206.6954   0.979
## bacArcicella sp.                   2258.5016  2364.4514  212.8589   0.955
## bacMicrobacterium \noxydans        1791.5725  2333.9711  214.0507   0.768
## bacPseudomonas \nprotogens         6424.0758  2344.7823  213.4976   2.740
## bacAll 10 bacteria                 7874.6934  2406.2248  214.3893   3.273
##                                    Pr(>|t|)    
## (Intercept)                         0.25865    
## pix0                                < 2e-16 ***
## bacFlavobacterium \nsuccinicans 1   0.89548    
## bacBosea \nmassiliensis             0.49368    
## bacAeromonas \nsalmonicida          0.93809    
## bacOhtaekwangia \nkoreensis         0.81138    
## bacFlavobacterium \nsuccinicans 2   0.44243    
## bacFalsiroseomonas \nstagni         0.75965    
## bacParasediminibacterium \npaludis  0.32882    
## bacArcicella sp.                    0.34056    
## bacMicrobacterium \noxydans         0.44357    
## bacPseudomonas \nprotogens          0.00667 ** 
## bacAll 10 bacteria                  0.00124 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 13 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```r
plot(lmer.grt.C2020)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-7.png)<!-- -->

```r
lmer.grt.C2020.2<-lmer(pix1 ~ pix0*bac + (1|edge) + (1|plate), data=df.C.2020)
anova(lmer.grt.C2020.2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##              Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)    
## pix0     3018862173 3018862173     1 207.78 60.0436 4.046e-13 ***
## bac       381160189   34650926    11 198.54  0.6892     0.748    
## pix0:bac  506580700   46052791    11 199.71  0.9160     0.526    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.C2020.2)
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + (1 | edge) + (1 | plate) + pix0:bac
##             npar  logLik    AIC    LRT Df Pr(>Chisq)    
## <none>        27 -2321.7 4697.3                         
## (1 | edge)    26 -2327.6 4707.2 11.871  1  0.0005703 ***
## (1 | plate)   26 -2325.3 4702.6  7.253  1  0.0070783 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.C2020.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 * bac + (1 | edge) + (1 | plate)
##    Data: df.C.2020
## 
## REML criterion at convergence: 4643.3
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.72877 -0.57352  0.00275  0.62103  3.11461 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept)  8715210 2952    
##  edge     (Intercept)  8488540 2914    
##  Residual             50277801 7091    
## Number of obs: 237, groups:  plate, 39; edge, 2
## 
## Fixed effects:
##                                           Estimate Std. Error         df
## (Intercept)                              6007.4638  6079.5059    53.3055
## pix0                                        4.1587     1.8724   208.0507
## bacFlavobacterium \nsuccinicans 1       -2933.8063  7886.4068   199.4464
## bacBosea \nmassiliensis                   536.5949  7140.6136   204.3446
## bacAeromonas \nsalmonicida                424.9411  8507.9317   207.1571
## bacOhtaekwangia \nkoreensis             -7654.1393  7688.4850   195.1687
## bacFlavobacterium \nsuccinicans 2       -1719.0370  8126.6166   197.4455
## bacFalsiroseomonas \nstagni              1721.8989  7421.0979   196.5983
## bacParasediminibacterium \npaludis      -1187.3873  6982.5276   202.9072
## bacArcicella sp.                         -628.1123  7901.5167   200.5321
## bacMicrobacterium \noxydans              8476.3254  7760.4207   201.0443
## bacPseudomonas \nprotogens              -6463.4287  7639.3435   188.8941
## bacAll 10 bacteria                       9294.7259 11476.4954   204.8952
## pix0:bacFlavobacterium \nsuccinicans 1      0.9321     2.9411   202.5169
## pix0:bacBosea \nmassiliensis                0.3381     2.3744   206.8263
## pix0:bacAeromonas \nsalmonicida            -0.2128     2.8181   209.9028
## pix0:bacOhtaekwangia \nkoreensis            2.3865     2.4691   199.1490
## pix0:bacFlavobacterium \nsuccinicans 2      1.2764     2.8473   195.9931
## pix0:bacFalsiroseomonas \nstagni           -0.3385     2.4110   200.4209
## pix0:bacParasediminibacterium \npaludis     1.2115     2.3126   201.0949
## pix0:bacArcicella sp.                       0.9628     2.5899   203.7609
## pix0:bacMicrobacterium \noxydans           -2.1995     2.5134   200.8961
## pix0:bacPseudomonas \nprotogens             4.6246     2.5878   193.5911
## pix0:bacAll 10 bacteria                    -0.8375     4.7039   204.8884
##                                         t value Pr(>|t|)  
## (Intercept)                               0.988   0.3275  
## pix0                                      2.221   0.0274 *
## bacFlavobacterium \nsuccinicans 1        -0.372   0.7103  
## bacBosea \nmassiliensis                   0.075   0.9402  
## bacAeromonas \nsalmonicida                0.050   0.9602  
## bacOhtaekwangia \nkoreensis              -0.996   0.3207  
## bacFlavobacterium \nsuccinicans 2        -0.212   0.8327  
## bacFalsiroseomonas \nstagni               0.232   0.8168  
## bacParasediminibacterium \npaludis       -0.170   0.8651  
## bacArcicella sp.                         -0.079   0.9367  
## bacMicrobacterium \noxydans               1.092   0.2760  
## bacPseudomonas \nprotogens               -0.846   0.3986  
## bacAll 10 bacteria                        0.810   0.4189  
## pix0:bacFlavobacterium \nsuccinicans 1    0.317   0.7516  
## pix0:bacBosea \nmassiliensis              0.142   0.8869  
## pix0:bacAeromonas \nsalmonicida          -0.076   0.9399  
## pix0:bacOhtaekwangia \nkoreensis          0.967   0.3350  
## pix0:bacFlavobacterium \nsuccinicans 2    0.448   0.6545  
## pix0:bacFalsiroseomonas \nstagni         -0.140   0.8885  
## pix0:bacParasediminibacterium \npaludis   0.524   0.6010  
## pix0:bacArcicella sp.                     0.372   0.7104  
## pix0:bacMicrobacterium \noxydans         -0.875   0.3826  
## pix0:bacPseudomonas \nprotogens           1.787   0.0755 .
## pix0:bacAll 10 bacteria                  -0.178   0.8589  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 24 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```r
plot(lmer.grt.C2020.2)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-8.png)<!-- -->

```r
#Wellspring 2020

lmer.grt.W2020<-lmer(pix1 ~ pix0 + bac + (1|edge) + (1|plate), data=df.W.2020)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.W2020)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##          Sum Sq    Mean Sq NumDF  DenDF  F value Pr(>F)    
## pix0 1.4664e+10 1.4664e+10     1 223.28 125.3584 <2e-16 ***
## bac  1.0798e+09 9.8160e+07    11 218.84   0.8392 0.6011    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.W2020)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + (1 | edge) + (1 | plate)
##             npar  logLik    AIC    LRT Df Pr(>Chisq)  
## <none>        16 -2444.1 4920.2                       
## (1 | edge)    15 -2445.9 4921.8 3.5107  1    0.06097 .
## (1 | plate)   15 -2445.2 4920.3 2.1062  1    0.14671  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.W2020)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 + bac + (1 | edge) + (1 | plate)
##    Data: df.W.2020
## 
## REML criterion at convergence: 4888.2
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.37937 -0.57901  0.08272  0.58502  2.42911 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev.
##  plate    (Intercept)   8054169  2838   
##  edge     (Intercept)   6735032  2595   
##  Residual             116974114 10815   
## Number of obs: 238, groups:  plate, 39; edge, 2
## 
## Fixed effects:
##                                    Estimate Std. Error         df t value
## (Intercept)                       5044.0552  3623.6227    10.5412   1.392
## pix0                                 7.2913     0.6512   223.2800  11.196
## bacSphingomonas \npituitosa 1     1484.5091  3560.5816   220.9672   0.417
## bacFlaviflagellibacter \ndeserti -1077.4527  3553.0940   221.5627  -0.303
## bacRhizobium \nrosettiformans      659.1365  3497.7890   219.4910   0.188
## bacRhizorhabdus \nwittichii 1     -581.6311  3490.8712   217.3165  -0.167
## bacRhizobium \ncapsici 1          1307.6764  3518.4586   223.0462   0.372
## bacPseudomonas \nprotogens        4856.7650  3527.2609   223.3201   1.377
## bacRhizorhabdus \nwittichii 2     3747.2298  3502.0644   220.2646   1.070
## bacRhizobium \ncapsici 2          4053.7219  3500.9568   219.8402   1.158
## bacSphingomonas \npituitosa 2     5636.4924  3502.5996   220.2029   1.609
## bacFervidobacterium \nriparium    1449.3299  3580.4124   224.0491   0.405
## bacAll 10 bacteria                4497.5803  3521.5576   222.8469   1.277
##                                  Pr(>|t|)    
## (Intercept)                         0.193    
## pix0                               <2e-16 ***
## bacSphingomonas \npituitosa 1       0.677    
## bacFlaviflagellibacter \ndeserti    0.762    
## bacRhizobium \nrosettiformans       0.851    
## bacRhizorhabdus \nwittichii 1       0.868    
## bacRhizobium \ncapsici 1            0.710    
## bacPseudomonas \nprotogens          0.170    
## bacRhizorhabdus \nwittichii 2       0.286    
## bacRhizobium \ncapsici 2            0.248    
## bacSphingomonas \npituitosa 2       0.109    
## bacFervidobacterium \nriparium      0.686    
## bacAll 10 bacteria                  0.203    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 13 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.W2020)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-9.png)<!-- -->

```r
lmer.grt.W2020.2<-lmer(pix1 ~ pix0*bac + (1|edge) + (1|plate), data=df.W.2020)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
anova(lmer.grt.W2020.2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##              Sum Sq    Mean Sq NumDF  DenDF  F value  Pr(>F)    
## pix0     1.4775e+10 1.4775e+10     1 211.70 131.9284 < 2e-16 ***
## bac      1.6500e+09 1.5000e+08    11 208.64   1.3394 0.20452    
## pix0:bac 2.3613e+09 2.1466e+08    11 209.37   1.9168 0.03875 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.grt.W2020.2)
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## pix1 ~ pix0 + bac + (1 | edge) + (1 | plate) + pix0:bac
##             npar  logLik    AIC    LRT Df Pr(>Chisq)  
## <none>        27 -2413.4 4880.9                       
## (1 | edge)    26 -2415.1 4882.2 3.2833  1    0.06999 .
## (1 | plate)   26 -2414.4 4880.9 1.9830  1    0.15908  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lmer.grt.W2020.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: pix1 ~ pix0 * bac + (1 | edge) + (1 | plate)
##    Data: df.W.2020
## 
## REML criterion at convergence: 4826.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.28221 -0.56032  0.09423  0.56351  2.38134 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev.
##  plate    (Intercept)   7642234  2764   
##  edge     (Intercept)   6663306  2581   
##  Residual             111989570 10583   
## Number of obs: 238, groups:  plate, 39; edge, 2
## 
## Fixed effects:
##                                         Estimate Std. Error         df t value
## (Intercept)                            16373.307   6417.353     74.104   2.551
## pix0                                       3.264      1.988    210.503   1.642
## bacSphingomonas \npituitosa 1         -10579.767   9420.055    208.015  -1.123
## bacFlaviflagellibacter \ndeserti      -12431.861   8732.440    209.582  -1.424
## bacRhizobium \nrosettiformans          -9259.457   9952.225    206.926  -0.930
## bacRhizorhabdus \nwittichii 1         -13969.398   8834.367    210.665  -1.581
## bacRhizobium \ncapsici 1               -6226.556   8943.044    210.893  -0.696
## bacPseudomonas \nprotogens            -24558.051   9927.076    209.231  -2.474
## bacRhizorhabdus \nwittichii 2         -12006.820   8930.869    207.534  -1.344
## bacRhizobium \ncapsici 2               -6135.807   9296.566    206.904  -0.660
## bacSphingomonas \npituitosa 2          -8303.000   8887.758    206.663  -0.934
## bacFervidobacterium \nriparium         10154.640  10639.530    211.572   0.954
## bacAll 10 bacteria                    -20119.553   9316.568    209.903  -2.160
## pix0:bacSphingomonas \npituitosa 1         4.354      3.275    208.873   1.330
## pix0:bacFlaviflagellibacter \ndeserti      4.051      2.759    212.134   1.468
## pix0:bacRhizobium \nrosettiformans         3.557      3.359    209.602   1.059
## pix0:bacRhizorhabdus \nwittichii 1         4.664      2.777    211.672   1.679
## pix0:bacRhizobium \ncapsici 1              2.675      2.898    209.965   0.923
## pix0:bacPseudomonas \nprotogens           11.388      3.625    208.802   3.142
## pix0:bacRhizorhabdus \nwittichii 2         5.690      2.990    208.654   1.903
## pix0:bacRhizobium \ncapsici 2              3.631      2.990    208.835   1.214
## pix0:bacSphingomonas \npituitosa 2         4.981      2.883    208.342   1.728
## pix0:bacFervidobacterium \nriparium       -1.299      2.981    211.650  -0.436
## pix0:bacAll 10 bacteria                    8.283      2.894    210.452   2.862
##                                       Pr(>|t|)   
## (Intercept)                            0.01279 * 
## pix0                                   0.10203   
## bacSphingomonas \npituitosa 1          0.26268   
## bacFlaviflagellibacter \ndeserti       0.15604   
## bacRhizobium \nrosettiformans          0.35325   
## bacRhizorhabdus \nwittichii 1          0.11532   
## bacRhizobium \ncapsici 1               0.48704   
## bacPseudomonas \nprotogens             0.01416 * 
## bacRhizorhabdus \nwittichii 2          0.18028   
## bacRhizobium \ncapsici 2               0.50998   
## bacSphingomonas \npituitosa 2          0.35129   
## bacFervidobacterium \nriparium         0.34096   
## bacAll 10 bacteria                     0.03194 * 
## pix0:bacSphingomonas \npituitosa 1     0.18509   
## pix0:bacFlaviflagellibacter \ndeserti  0.14362   
## pix0:bacRhizobium \nrosettiformans     0.29086   
## pix0:bacRhizorhabdus \nwittichii 1     0.09461 . 
## pix0:bacRhizobium \ncapsici 1          0.35706   
## pix0:bacPseudomonas \nprotogens        0.00192 **
## pix0:bacRhizorhabdus \nwittichii 2     0.05843 . 
## pix0:bacRhizobium \ncapsici 2          0.22607   
## pix0:bacSphingomonas \npituitosa 2     0.08554 . 
## pix0:bacFervidobacterium \nriparium    0.66337   
## pix0:bacAll 10 bacteria                0.00463 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 24 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
plot(lmer.grt.W2020.2)
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-10.png)<!-- -->

```r
##########################################################################################

# Onto figure generation - plotting Log Response Ratios of bacterial effects on plant growth

# Churchill 2019

# Reference for calculating LRR's: Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis of response ratios in experimental ecology. Ecology, 80(4), 1150-1156. Equation 1

# Let's do this more straightforwardly. For each bacterial treatment, we need count, mean, and SD.

melt.C.2019<-melt(df.C.2019, id.vars=c("bac"), measure.vars= "Grt", na.rm = T)
Sum.C.2019<- ddply(melt.C.2019, c("bac","variable"), dplyr::summarise, mean = mean(value), sd = sd(value), count=n(),  sem = sd(value)/sqrt(length(value)))

lrr.C.2019<-Sum.C.2019[2:12,]
lrr.C.2019$LRR<-log(lrr.C.2019$mean) - log(Sum.C.2019$mean[1])
lrr.C.2019$var<-lrr.C.2019$sd^2/(lrr.C.2019$count*lrr.C.2019$mean^2) + Sum.C.2019$sd[1]^2/(Sum.C.2019$count[1]*Sum.C.2019$mean[1]^2)
lrr.C.2019$CI<-qnorm(0.975)*sqrt(lrr.C.2019$var)
lrr.C.2019$SE<-sqrt(lrr.C.2019$var)/sqrt(lrr.C.2019$count)

# Kelso 2019

melt.K.2019<-melt(df.K.2019, id.vars=c("bac"), measure.vars= "Grt", na.rm = T)
Sum.K.2019<- ddply(melt.K.2019, c("bac","variable"), dplyr::summarise, mean = mean(value), sd = sd(value), count=n(), sem = sd(value)/sqrt(length(value)))

lrr.K.2019<-Sum.K.2019[2:12,]
lrr.K.2019$LRR<-log(lrr.K.2019$mean) - log(Sum.K.2019$mean[1])
lrr.K.2019$var<-lrr.K.2019$sd^2/(lrr.K.2019$count*lrr.K.2019$mean^2) + Sum.K.2019$sd[1]^2/(Sum.K.2019$count[1]*Sum.K.2019$mean[1]^2)
lrr.K.2019$CI<-qnorm(0.975)*sqrt(lrr.K.2019$var)
lrr.K.2019$SE<-sqrt(lrr.K.2019$var)/sqrt(lrr.K.2019$count)
# Wellspring 2019

melt.W.2019<-melt(df.W.2019, id.vars=c("bac"), measure.vars= "Grt", na.rm = T)
Sum.W.2019<- ddply(melt.W.2019, c("bac","variable"), dplyr::summarise, mean = mean(value), sd = sd(value), count=n(), sem = sd(value)/sqrt(length(value)))

lrr.W.2019<-Sum.W.2019[2:12,]
lrr.W.2019$LRR<-log(lrr.W.2019$mean) - log(Sum.W.2019$mean[1])
lrr.W.2019$var<-lrr.W.2019$sd^2/(lrr.W.2019$count*lrr.W.2019$mean^2) + Sum.W.2019$sd[1]^2/(Sum.W.2019$count[1]*Sum.W.2019$mean[1]^2)
lrr.W.2019$CI<-qnorm(0.975)*sqrt(lrr.W.2019$var)
lrr.W.2019$SE<-sqrt(lrr.W.2019$var)/sqrt(lrr.W.2019$count)

# Churchill 2020

melt.C.2020<-melt(df.C.2020, id.vars=c("bac"), measure.vars= "Grt", na.rm = T)
Sum.C.2020<- ddply(melt.C.2020, c("bac","variable"), dplyr::summarise, mean = mean(value), sd = sd(value), count=n(), sem = sd(value)/sqrt(length(value)))

lrr.C.2020<-Sum.C.2020[2:12,]
lrr.C.2020$LRR<-log(lrr.C.2020$mean) - log(Sum.C.2020$mean[1])
lrr.C.2020$var<-lrr.C.2020$sd^2/(lrr.C.2020$count*lrr.C.2020$mean^2) + Sum.C.2020$sd[1]^2/(Sum.C.2020$count[1]*Sum.C.2020$mean[1]^2)
lrr.C.2020$CI<-qnorm(0.975)*sqrt(lrr.C.2020$var)
lrr.C.2020$SE<-sqrt(lrr.C.2020$var)/sqrt(lrr.C.2020$count)

#Wellspring
melt.W.2020<-melt(df.W.2020, id.vars=c("bac"), measure.vars= "Grt", na.rm = T)
Sum.W.2020<- ddply(melt.W.2020, c("bac","variable"), dplyr::summarise, mean = mean(value), sd = sd(value), count=n(), sem = sd(value)/sqrt(length(value)))

lrr.W.2020<-Sum.W.2020[2:12,]
lrr.W.2020$LRR<-log(lrr.W.2020$mean) - log(Sum.W.2020$mean[1])
lrr.W.2020$var<-lrr.W.2020$sd^2/(lrr.W.2020$count*lrr.W.2020$mean^2) + Sum.W.2020$sd[1]^2/(Sum.W.2020$count[1]*Sum.W.2020$mean[1]^2)
lrr.W.2020$CI<-qnorm(0.975)*sqrt(lrr.W.2020$var)
lrr.W.2020$SE<-sqrt(lrr.W.2020$var)/sqrt(lrr.W.2020$count)

#OK, let's get to plotting!

pred.C.2019<-mean(lrr.C.2019$LRR[1:10])

plot1.C1 <- ggplot(lrr.C.2019,aes(y=LRR,x=bac,colour=bac))+geom_errorbar(aes(ymin=LRR-SE,ymax=LRR+SE),width=0.5)+geom_point(size=3)
plot1.C1 <- plot1.C1 + geom_hline(yintercept=pred.C.2019, linetype="dashed", color = "blue3",size=1)
plot1.C1 <- plot1.C1 + theme_classic() + ggtitle("Churchill, 2019") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot1.C1 <- plot1.C1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot1.C1 <- plot1.C1 + theme(legend.position = "none")
plot1.C1 <- plot1.C1 + labs(x= element_blank()) 
plot1.C1 <- plot1.C1 + labs(y= element_blank())
plot1.C1 <- plot1.C1 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))
plot1.C1 <- plot1.C1 + geom_hline(yintercept=0, linetype="dashed", color = "red3",size=1)

#Kelso 2019

pred.K.2019<-mean(lrr.K.2019$LRR[1:10])

plot1.K1 <- ggplot(lrr.K.2019,aes(y=LRR,x=bac,colour=bac))+geom_errorbar(aes(ymin=LRR-SE,ymax=LRR+SE),width=0.5)+geom_point(size=3)
plot1.K1 <- plot1.K1 + geom_hline(yintercept=pred.K.2019, linetype="dashed", color = "blue3",size=1)
plot1.K1 <- plot1.K1 + theme_classic() + ggtitle("Kelso, 2019") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot1.K1 <- plot1.K1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot1.K1 <- plot1.K1 + theme(legend.position = "none")
plot1.K1 <- plot1.K1 + labs(x= element_blank()) 
plot1.K1 <- plot1.K1 + labs(y= element_blank())
plot1.K1 <- plot1.K1 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))
plot1.K1 <- plot1.K1 + geom_hline(yintercept=0, linetype="dashed", color = "red3",size=1)

#Wellspring 2019

pred.W.2019<-mean(lrr.W.2019$LRR[1:10])

plot1.W1 <- ggplot(lrr.W.2019,aes(y=LRR,x=bac,colour=bac))+geom_errorbar(aes(ymin=LRR-SE,ymax=LRR+SE),width=0.5)+geom_point(size=3)
plot1.W1 <- plot1.W1 + geom_hline(yintercept=pred.W.2019, linetype="dashed", color = "blue3",size=1)
plot1.W1 <- plot1.W1 + theme_classic() + ggtitle("Wellspring, 2019") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot1.W1 <- plot1.W1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot1.W1 <- plot1.W1 + theme(legend.position = "none")
plot1.W1 <- plot1.W1 + labs(x= element_blank()) 
plot1.W1 <- plot1.W1 + labs(y= element_blank())
plot1.W1 <- plot1.W1 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))
plot1.W1 <- plot1.W1 + geom_hline(yintercept=0, linetype="dashed", color = "red3",size=1)

# Churchill 2020

pred.C.2020<-mean(lrr.C.2020$LRR[1:10])

plot1.C2 <- ggplot(lrr.C.2020,aes(y=LRR,x=bac,colour=bac))+geom_errorbar(aes(ymin=LRR-SE,ymax=LRR+SE),width=0.5)+geom_point(size=3)
plot1.C2 <- plot1.C2 + geom_hline(yintercept=pred.C.2020, linetype="dashed", color = "blue3",size=1)
plot1.C2 <- plot1.C2 + theme_classic() + ggtitle("Churchill, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot1.C2 <- plot1.C2 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot1.C2 <- plot1.C2 + theme(legend.position = "none")
plot1.C2 <- plot1.C2 + labs(x= element_blank()) 
plot1.C2 <- plot1.C2 + labs(y= element_blank())
plot1.C2 <- plot1.C2 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))
plot1.C2 <- plot1.C2 + geom_hline(yintercept=0, linetype="dashed", color = "red3",size=1)

# Wellspring 2020

pred.W.2020<-mean(lrr.W.2020$LRR[1:10])

plot1.W2 <- ggplot(lrr.W.2020,aes(y=LRR,x=bac,colour=bac))+geom_errorbar(aes(ymin=LRR-SE,ymax=LRR+SE),width=0.5)+geom_point(size=3)
plot1.W2 <- plot1.W2 + geom_hline(yintercept=pred.W.2020, linetype="dashed", color = "blue3",size=1)
plot1.W2 <- plot1.W2 + theme_classic() + ggtitle("Wellspring, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot1.W2 <- plot1.W2 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot1.W2 <- plot1.W2 + theme(legend.position = "none")
plot1.W2 <- plot1.W2 + labs(x= element_blank()) 
plot1.W2 <- plot1.W2 + labs(y= element_blank())
plot1.W2 <- plot1.W2 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))
plot1.W2 <- plot1.W2 + geom_hline(yintercept=0, linetype="dashed", color = "red3",size=1)

#Legend plot
#Create dead column for individual bacteria
lrr.C.2019$legend.bac<-c(rep('Individual bacteria',10),'All 10 bacteria')
lrr.C.2019$legend.bac<-as.factor(lrr.C.2019$legend.bac)
factor(lrr.C.2019$legend.bac, levels = c("Individual bacteria","All 10 bacteria"))
```

```
##  [1] Individual bacteria Individual bacteria Individual bacteria
##  [4] Individual bacteria Individual bacteria Individual bacteria
##  [7] Individual bacteria Individual bacteria Individual bacteria
## [10] Individual bacteria All 10 bacteria    
## Levels: Individual bacteria All 10 bacteria
```

```r
lrr.C.2019$legend.line<-c(rep('No effect of inoculation',6),rep('Mean effect of individual bacteria',5))
lrr.C.2019$legend.line<-as.factor(lrr.C.2019$legend.line)

leg.plot <- ggplot(lrr.C.2019, aes(x = LRR, y = var, colour = legend.bac)) + geom_point(size=3)
leg.plot <- leg.plot + scale_colour_manual(name="Bacterial inocula",values=c("blue3","#000000"))
leg.plot <- leg.plot + theme_classic()
leg.plot <- leg.plot + theme(
  legend.title=element_text(size=14, face='bold'),
  legend.text=element_text(size=12, face='bold')
)

leg.1 <- get_legend(leg.plot)  

leg.plot.2 <- ggplot(lrr.C.2019, aes(x = LRR, y = var, colour = legend.line)) + geom_point(size=3) + geom_line(linetype="dashed",size=1)
leg.plot.2 <- leg.plot.2 + scale_colour_manual(name="Line",values=c("blue3","red3"))
leg.plot.2 <- leg.plot.2 + theme_classic() 
leg.plot.2 <- leg.plot.2 + theme(
  legend.title=element_text(size=14, face='bold'),
  legend.text=element_text(size=12, face='bold')
)

leg.2 <- get_legend(leg.plot.2)  

plot.leg<-plot_grid(leg.1, leg.2, ncol=1,nrow=2, align='v')

# Generate figure 1

x.grob_bac <- textGrob(expression(bold("Bacterial strain")), gp=gpar(fontsize=12))
y.grob_lrrdw <- textGrob(expression(bold("                     Duckweed growth (LRR)")), gp=gpar(fontsize=12), rot=90)

# Re-scale y axes
plot1.C1 <- plot1.C1 + ylim(-2.5,1)
plot1.K1 <- plot1.K1 + ylim(-2.5,1)
plot1.W1 <- plot1.W1 + ylim(-2.5,1)

plot1.C2 <- plot1.C2 + ylim(-0.2,0.4)
plot1.W2 <- plot1.W2 + ylim(-0.2,0.4)

plot1.2019<-plot_grid(plot1.C1, plot1.K1, plot1.W1, ncol=3,nrow=1, align='h', labels=c('A','B','C'))
plot1.2020<-plot_grid(plot1.C2, plot1.W2, plot.leg, ncol=3,nrow=1, align='h', labels=c('D','E'))

plot1<-plot_grid(plot1.2019,plot1.2020, ncol=1, nrow=2, align='v')
grid.arrange(arrangeGrob(plot1, bottom = x.grob_bac, left = y.grob_lrrdw))
```

![](Ch1_final_fuller_code_files/figure-html/Duckweed growth-11.png)<!-- -->

```r
# PDF , 10 x 12.5 " works well
```



SECTION 4: Microbial growth models and figures


```r
#Churchill 2020

# Model 1 : we will exclude control plants, and use All10_pltN as our reference treatment.
df.C.2020$bac2<-df.C.2020$bac
df.C.2020 <- within(df.C.2020, bac2 <- relevel(bac2, ref = 12))

lmer.abs.C<-lmer(logabs~bac2*plt + (1|edge) + (1|plate), data=df.C.2020)
```

```
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```
## boundary (singular) fit: see help('isSingular')
```

```r
anova(lmer.abs.C)
```

```
## Missing cells for: bac2Control:pltN.  
## Interpret type III hypotheses with care.
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
## bac2     182.795  16.618    11 207.94  18.657 < 2.2e-16 ***
## plt      168.266 168.266     1 211.50 188.910 < 2.2e-16 ***
## bac2:plt  42.318   4.232    10 206.19   4.751 3.938e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.abs.C)
```

```
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## logabs ~ bac2 + plt + (1 | edge) + (1 | plate) + bac2:plt
##             npar  logLik    AIC    LRT Df Pr(>Chisq)  
## <none>        26 -321.44 694.89                       
## (1 | edge)    25 -321.44 692.89 0.0000  1    1.00000  
## (1 | plate)   25 -324.30 698.60 5.7138  1    0.01683 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
plot(lmer.abs.C)
```

![](Ch1_final_fuller_code_files/figure-html/Microbial growth models-1.png)<!-- -->

```r
summary(lmer.abs.C)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: logabs ~ bac2 * plt + (1 | edge) + (1 | plate)
##    Data: df.C.2020
## 
## REML criterion at convergence: 642.9
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.0544 -0.4200  0.0686  0.5012  3.3798 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept) 0.08433  0.2904  
##  edge     (Intercept) 0.00000  0.0000  
##  Residual             0.89072  0.9438  
## Number of obs: 235, groups:  plate, 19; edge, 2
## 
## Fixed effects:
##                                          Estimate Std. Error       df t value
## (Intercept)                                7.8529     0.3119 207.3425  25.180
## bac2Control                               -1.1012     0.4505 209.4154  -2.444
## bac2Flavobacterium \nsuccinicans 1        -3.1387     0.4250 210.4944  -7.385
## bac2Bosea \nmassiliensis                  -3.1418     0.4581 208.1761  -6.858
## bac2Aeromonas \nsalmonicida               -1.5574     0.4164 210.3544  -3.741
## bac2Ohtaekwangia \nkoreensis              -2.5924     0.4341 209.6444  -5.972
## bac2Flavobacterium \nsuccinicans 2        -2.8365     0.4127 207.1377  -6.873
## bac2Falsiroseomonas \nstagni              -3.2959     0.4003 208.2845  -8.235
## bac2Parasediminibacterium \npaludis       -0.1782     0.4979 207.4435  -0.358
## bac2Arcicella sp.                         -2.2163     0.4141 208.4927  -5.352
## bac2Microbacterium \noxydans              -3.4730     0.4152 209.3080  -8.364
## bac2Pseudomonas \nprotogens                0.1757     0.4419 206.4491   0.398
## pltY                                       0.6347     0.4534 203.1265   1.400
## bac2Flavobacterium \nsuccinicans 1:pltY    2.0591     0.6518 206.8515   3.159
## bac2Bosea \nmassiliensis:pltY              1.7793     0.6365 208.3730   2.795
## bac2Aeromonas \nsalmonicida:pltY           0.4964     0.5865 204.0302   0.846
## bac2Ohtaekwangia \nkoreensis:pltY          0.4037     0.6481 208.0255   0.623
## bac2Flavobacterium \nsuccinicans 2:pltY    2.1655     0.6197 203.9089   3.494
## bac2Falsiroseomonas \nstagni:pltY          2.3367     0.5877 206.0861   3.976
## bac2Parasediminibacterium \npaludis:pltY  -0.1190     0.6884 204.2244  -0.173
## bac2Arcicella sp.:pltY                     1.5378     0.6050 204.0039   2.542
## bac2Microbacterium \noxydans:pltY          2.1594     0.6188 202.9168   3.490
## bac2Pseudomonas \nprotogens:pltY           0.5143     0.6521 205.8084   0.789
##                                          Pr(>|t|)    
## (Intercept)                               < 2e-16 ***
## bac2Control                              0.015335 *  
## bac2Flavobacterium \nsuccinicans 1       3.51e-12 ***
## bac2Bosea \nmassiliensis                 7.79e-11 ***
## bac2Aeromonas \nsalmonicida              0.000237 ***
## bac2Ohtaekwangia \nkoreensis             9.89e-09 ***
## bac2Flavobacterium \nsuccinicans 2       7.24e-11 ***
## bac2Falsiroseomonas \nstagni             1.97e-14 ***
## bac2Parasediminibacterium \npaludis      0.720820    
## bac2Arcicella sp.                        2.29e-07 ***
## bac2Microbacterium \noxydans             8.51e-15 ***
## bac2Pseudomonas \nprotogens              0.691288    
## pltY                                     0.163076    
## bac2Flavobacterium \nsuccinicans 1:pltY  0.001819 ** 
## bac2Bosea \nmassiliensis:pltY            0.005667 ** 
## bac2Aeromonas \nsalmonicida:pltY         0.398293    
## bac2Ohtaekwangia \nkoreensis:pltY        0.534020    
## bac2Flavobacterium \nsuccinicans 2:pltY  0.000583 ***
## bac2Falsiroseomonas \nstagni:pltY        9.70e-05 ***
## bac2Parasediminibacterium \npaludis:pltY 0.862932    
## bac2Arcicella sp.:pltY                   0.011777 *  
## bac2Microbacterium \noxydans:pltY        0.000593 ***
## bac2Pseudomonas \nprotogens:pltY         0.431223    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 23 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
#Model 2 : we will exclude wells without plants to estimate variation in bacterial cell density when co-cultured with plants.
df.C.2020.plt<-subset(df.C.2020,df.C.2020$plt=="Y")
lmer.abs.C.2<-lmer(logabs~bac2 + (1|edge) + (1|plate), data=df.C.2020.plt)
```

```
## boundary (singular) fit: see help('isSingular')
```

```r
summary(lmer.abs.C.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: logabs ~ bac2 + (1 | edge) + (1 | plate)
##    Data: df.C.2020.plt
## 
## REML criterion at convergence: 300.5
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.3001 -0.2930  0.1813  0.4878  1.7692 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept) 0.000000 0.00000 
##  edge     (Intercept) 0.005539 0.07443 
##  Residual             0.750792 0.86648 
## Number of obs: 119, groups:  plate, 19; edge, 2
## 
## Fixed effects:
##                                     Estimate Std. Error       df t value
## (Intercept)                           8.4247     0.3127  36.0592  26.943
## bac2Control                          -1.0829     0.4028 106.5647  -2.688
## bac2Flavobacterium \nsuccinicans 1   -0.9925     0.4485 106.0080  -2.213
## bac2Bosea \nmassiliensis             -1.2532     0.3964 106.8569  -3.162
## bac2Aeromonas \nsalmonicida          -0.9646     0.3798 106.9443  -2.540
## bac2Ohtaekwangia \nkoreensis         -2.0882     0.4338 106.9756  -4.813
## bac2Flavobacterium \nsuccinicans 2   -0.5413     0.4215 106.8555  -1.284
## bac2Falsiroseomonas \nstagni         -0.9094     0.3897 106.8039  -2.333
## bac2Parasediminibacterium \npaludis  -0.2400     0.4338 106.9756  -0.553
## bac2Arcicella sp.                    -0.6406     0.4026 106.0358  -1.591
## bac2Microbacterium \noxydans         -1.3277     0.4215 106.8555  -3.150
## bac2Pseudomonas \nprotogens           0.6764     0.4338 106.9756   1.559
##                                     Pr(>|t|)    
## (Intercept)                          < 2e-16 ***
## bac2Control                          0.00834 ** 
## bac2Flavobacterium \nsuccinicans 1   0.02902 *  
## bac2Bosea \nmassiliensis             0.00204 ** 
## bac2Aeromonas \nsalmonicida          0.01253 *  
## bac2Ohtaekwangia \nkoreensis         4.9e-06 ***
## bac2Flavobacterium \nsuccinicans 2   0.20176    
## bac2Falsiroseomonas \nstagni         0.02150 *  
## bac2Parasediminibacterium \npaludis  0.58121    
## bac2Arcicella sp.                    0.11460    
## bac2Microbacterium \noxydans         0.00212 ** 
## bac2Pseudomonas \nprotogens          0.12194    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) bc2Cnt bc2Fs1 bc2Bsm bc2Ars bc2Ohk bc2Fs2 bc2Fls bc2Prp
## bac2Control -0.749                                                        
## bc2Flvbcts1 -0.670  0.520                                                 
## bc2Bsmsslns -0.764  0.590  0.528                                          
## bc2Armnsslm -0.796  0.615  0.551  0.627                                   
## bc2Ohtkwngk -0.698  0.539  0.483  0.549  0.572                            
## bc2Flvbcts2 -0.717  0.554  0.497  0.565  0.589  0.516                     
## bc2Flsrsmns -0.775  0.599  0.537  0.610  0.636  0.557  0.573              
## bc2Prsdmnbp -0.698  0.539  0.483  0.549  0.572  0.501  0.516  0.557       
## bc2Arcclsp. -0.745  0.578  0.520  0.588  0.613  0.537  0.553  0.598  0.537
## bc2Mcrbctro -0.717  0.554  0.497  0.565  0.589  0.516  0.530  0.573  0.516
## bc2Psdmnspr -0.698  0.539  0.483  0.549  0.572  0.501  0.516  0.557  0.501
##             bc2As. bc2Mco
## bac2Control              
## bc2Flvbcts1              
## bc2Bsmsslns              
## bc2Armnsslm              
## bc2Ohtkwngk              
## bc2Flvbcts2              
## bc2Flsrsmns              
## bc2Prsdmnbp              
## bc2Arcclsp.              
## bc2Mcrbctro  0.553       
## bc2Psdmnspr  0.537  0.516
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
#Wellspring 2020

df.W.2020$bac2<-df.W.2020$bac
df.W.2020 <- within(df.W.2020, bac2 <- relevel(bac2, ref = 12))

lmer.abs.W<-lmer(logabs~bac2*plt + (1|edge) + (1|plate), data=df.W.2020)
```

```
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```r
anova(lmer.abs.W, type=3)
```

```
## Missing cells for: bac2Control:pltN.  
## Interpret type III hypotheses with care.
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##           Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
## bac2     132.509  12.046    11 190.99  15.1603 < 2.2e-16 ***
## plt      220.099 220.099     1 193.00 276.9952 < 2.2e-16 ***
## bac2:plt  35.861   3.586    10 191.28   4.5131 9.792e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
ranova(lmer.abs.W)
```

```
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```
## ANOVA-like table for random-effects: Single term deletions
## 
## Model:
## logabs ~ bac2 + plt + (1 | edge) + (1 | plate) + bac2:plt
##             npar  logLik    AIC    LRT Df Pr(>Chisq)  
## <none>        26 -288.66 629.31                       
## (1 | edge)    25 -290.93 631.85 4.5377  1    0.03316 *
## (1 | plate)   25 -289.86 629.73 2.4112  1    0.12047  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
plot(lmer.abs.W)
```

![](Ch1_final_fuller_code_files/figure-html/Microbial growth models-2.png)<!-- -->

```r
summary(lmer.abs.W)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: logabs ~ bac2 * plt + (1 | edge) + (1 | plate)
##    Data: df.W.2020
## 
## REML criterion at convergence: 577.3
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.1033 -0.4839  0.0301  0.5717  3.4841 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept) 0.05782  0.2405  
##  edge     (Intercept) 0.06233  0.2497  
##  Residual             0.79460  0.8914  
## Number of obs: 220, groups:  plate, 19; edge, 2
## 
## Fixed effects:
##                                        Estimate Std. Error       df t value
## (Intercept)                              8.2385     0.3292   8.9595  25.024
## bac2Control                             -1.2316     0.4960 188.0840  -2.483
## bac2Sphingomonas \npituitosa 1          -1.7858     0.4141 196.0457  -4.312
## bac2Flaviflagellibacter \ndeserti       -4.4645     0.4123 194.6515 -10.828
## bac2Rhizobium \nrosettiformans          -2.8909     0.3951 186.4410  -7.318
## bac2Rhizorhabdus \nwittichii 1          -3.5581     0.3654 189.9517  -9.738
## bac2Rhizobium \ncapsici 1               -2.8390     0.4200 188.2331  -6.759
## bac2Pseudomonas \nprotogens             -1.0604     0.4395 190.3181  -2.413
## bac2Rhizorhabdus \nwittichii 2          -3.5998     0.4261 194.5314  -8.447
## bac2Rhizobium \ncapsici 2               -3.1476     0.4096 192.9921  -7.685
## bac2Sphingomonas \npituitosa 2          -2.4406     0.3957 188.9233  -6.168
## bac2Fervidobacterium \nriparium         -2.7748     0.3972 188.7894  -6.986
## pltY                                     0.5586     0.4917 188.4010   1.136
## bac2Sphingomonas \npituitosa 1:pltY      0.6294     0.6510 195.0919   0.967
## bac2Flaviflagellibacter \ndeserti:pltY   2.0737     0.6742 189.3303   3.076
## bac2Rhizobium \nrosettiformans:pltY      1.9954     0.6577 192.8808   3.034
## bac2Rhizorhabdus \nwittichii 1:pltY      2.0797     0.6222 188.4762   3.342
## bac2Rhizobium \ncapsici 1:pltY           0.7970     0.6317 184.3681   1.262
## bac2Pseudomonas \nprotogens:pltY         0.4657     0.6655 193.4255   0.700
## bac2Rhizorhabdus \nwittichii 2:pltY      2.2963     0.6707 192.6101   3.424
## bac2Rhizobium \ncapsici 2:pltY           2.6764     0.6315 187.7431   4.238
## bac2Sphingomonas \npituitosa 2:pltY      2.3314     0.6338 189.7456   3.678
## bac2Fervidobacterium \nriparium:pltY     2.2067     0.6318 190.7527   3.493
##                                        Pr(>|t|)    
## (Intercept)                            1.34e-09 ***
## bac2Control                            0.013904 *  
## bac2Sphingomonas \npituitosa 1         2.56e-05 ***
## bac2Flaviflagellibacter \ndeserti       < 2e-16 ***
## bac2Rhizobium \nrosettiformans         7.29e-12 ***
## bac2Rhizorhabdus \nwittichii 1          < 2e-16 ***
## bac2Rhizobium \ncapsici 1              1.69e-10 ***
## bac2Pseudomonas \nprotogens            0.016777 *  
## bac2Rhizorhabdus \nwittichii 2         6.87e-15 ***
## bac2Rhizobium \ncapsici 2              7.54e-13 ***
## bac2Sphingomonas \npituitosa 2         4.12e-09 ***
## bac2Fervidobacterium \nriparium        4.70e-11 ***
## pltY                                   0.257368    
## bac2Sphingomonas \npituitosa 1:pltY    0.334874    
## bac2Flaviflagellibacter \ndeserti:pltY 0.002411 ** 
## bac2Rhizobium \nrosettiformans:pltY    0.002748 ** 
## bac2Rhizorhabdus \nwittichii 1:pltY    0.001002 ** 
## bac2Rhizobium \ncapsici 1:pltY         0.208665    
## bac2Pseudomonas \nprotogens:pltY       0.484921    
## bac2Rhizorhabdus \nwittichii 2:pltY    0.000754 ***
## bac2Rhizobium \ncapsici 2:pltY         3.53e-05 ***
## bac2Sphingomonas \npituitosa 2:pltY    0.000306 ***
## bac2Fervidobacterium \nriparium:pltY   0.000593 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Correlation matrix not shown by default, as p = 23 > 12.
## Use print(x, correlation=TRUE)  or
##     vcov(x)        if you need it
```

```
## fit warnings:
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

```r
df.W.2020.plt<-subset(df.W.2020,df.W.2020$plt=="Y")
lmer.abs.W.2<-lmer(logabs~bac2 + (1|edge) + (1|plate), data=df.W.2020.plt)
```

```
## boundary (singular) fit: see help('isSingular')
```

```r
summary(lmer.abs.W.2)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: logabs ~ bac2 + (1 | edge) + (1 | plate)
##    Data: df.W.2020.plt
## 
## REML criterion at convergence: 304.3
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.1444 -0.3810  0.1059  0.5948  1.6410 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  plate    (Intercept) 0.00000  0.0000  
##  edge     (Intercept) 0.09692  0.3113  
##  Residual             0.85191  0.9230  
## Number of obs: 115, groups:  plate, 19; edge, 2
## 
## Fixed effects:
##                                   Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)                         8.8390     0.4764  13.3850  18.555 6.14e-11
## bac2Control                        -1.2660     0.5084 102.4079  -2.490  0.01438
## bac2Sphingomonas \npituitosa 1     -1.2266     0.5135 102.8634  -2.389  0.01873
## bac2Flaviflagellibacter \ndeserti  -2.4078     0.5459 102.6506  -4.411 2.55e-05
## bac2Rhizobium \nrosettiformans     -0.9623     0.5305 102.5527  -1.814  0.07259
## bac2Rhizorhabdus \nwittichii 1     -1.4986     0.5164 102.2299  -2.902  0.00454
## bac2Rhizobium \ncapsici 1          -1.9966     0.4875 102.2745  -4.096 8.43e-05
## bac2Pseudomonas \nprotogens        -0.6988     0.5045 102.7877  -1.385  0.16898
## bac2Rhizorhabdus \nwittichii 2     -1.3609     0.5267 102.0730  -2.584  0.01118
## bac2Rhizobium \ncapsici 2          -0.5299     0.4949 102.5108  -1.071  0.28680
## bac2Sphingomonas \npituitosa 2     -0.1717     0.5074 102.9467  -0.338  0.73577
## bac2Fervidobacterium \nriparium    -0.4298     0.4981 102.0435  -0.863  0.39020
##                                      
## (Intercept)                       ***
## bac2Control                       *  
## bac2Sphingomonas \npituitosa 1    *  
## bac2Flaviflagellibacter \ndeserti ***
## bac2Rhizobium \nrosettiformans    .  
## bac2Rhizorhabdus \nwittichii 1    ** 
## bac2Rhizobium \ncapsici 1         ***
## bac2Pseudomonas \nprotogens          
## bac2Rhizorhabdus \nwittichii 2    *  
## bac2Rhizobium \ncapsici 2            
## bac2Sphingomonas \npituitosa 2       
## bac2Fervidobacterium \nriparium      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) bc2Cnt bc2Sp1 bc2Fld bc2Rhr bc2Rw1 bc2Rc1 bc2Psp bc2Rw2
## bac2Control -0.724                                                        
## bc2Sphngmp1 -0.730  0.671                                                 
## bc2Flvflgld -0.682  0.629  0.633                                          
## bc2Rhzbmrst -0.698  0.645  0.648  0.606                                   
## bc2Rhzrhbw1 -0.707  0.657  0.656  0.615  0.632                            
## bc2Rhzbmcp1 -0.750  0.696  0.696  0.652  0.670  0.683                     
## bc2Psdmnspr -0.740  0.681  0.686  0.642  0.657  0.667  0.707              
## bc2Rhzrhbw2 -0.687  0.641  0.638  0.599  0.615  0.630  0.667  0.648       
## bc2Rhzbmcp2 -0.746  0.690  0.692  0.648  0.664  0.676  0.716  0.702  0.659
## bc2Sphngmp2 -0.741  0.681  0.688  0.642  0.658  0.665  0.705  0.697  0.646
## bc2Frvdbctr -0.724  0.676  0.672  0.631  0.649  0.665  0.705  0.683  0.651
##             bc2Rc2 bc2Sp2
## bac2Control              
## bc2Sphngmp1              
## bc2Flvflgld              
## bc2Rhzbmrst              
## bc2Rhzrhbw1              
## bc2Rhzbmcp1              
## bc2Psdmnspr              
## bc2Rhzrhbw2              
## bc2Rhzbmcp2              
## bc2Sphngmp2  0.702       
## bc2Frvdbctr  0.695  0.681
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
# Generate Figure 2

melt.C.abs<-melt(df.C.2020, id.vars=c("bac","plt"), measure.vars= "abs", na.rm = T)
Sum.C.abs<- ddply(melt.C.abs, c("bac","plt","variable"), summarise, mean = mean(value), sd = sd(value), count=n(), sem = sd(value)/sqrt(length(value)))

# Re-order by bacterial population size

Sum.C.abs$bac <- factor(Sum.C.abs$bac, levels = c("Control","Microbacterium \noxydans","Bosea \nmassiliensis","Flavobacterium \nsuccinicans 1","Falsiroseomonas \nstagni","Ohtaekwangia \nkoreensis","Flavobacterium \nsuccinicans 2","Aeromonas \nsalmonicida","Arcicella sp.","Parasediminibacterium \npaludis","Pseudomonas \nprotogens","All 10 bacteria"))
Sum.C.abs$CI<-qnorm(0.975)*Sum.C.abs$sd/sqrt(Sum.C.abs$count)

Sum.C.plt<-subset(Sum.C.abs, Sum.C.abs$plt=="Y")
Sum.C.plt<-Sum.C.plt[2:12,]
Sum.C.plt$CI<-qnorm(0.975)*Sum.C.plt$sd/sqrt(Sum.C.plt$count)

Sum.C.bac<-subset(Sum.C.abs, Sum.C.abs$plt=="N")
Sum.C.bac$CI<-qnorm(0.975)*Sum.C.bac$sd/sqrt(Sum.C.bac$count)

pred.C.plt<-mean(Sum.C.plt$mean[1:10])
pred.C.bac<-mean(Sum.C.bac$mean[1:10])

plot2.C1 <- ggplot(Sum.C.bac,aes(y=mean,x=bac,colour=bac))+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=0.5)+geom_point(size=3)
plot2.C1 <- plot2.C1 + geom_hline(yintercept=pred.C.bac, linetype="dashed", color = "blue3",size=1)
plot2.C1 <- plot2.C1 + theme_classic() + ggtitle("Churchill, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot2.C1 <- plot2.C1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot2.C1 <- plot2.C1 + theme(legend.position = "none")
plot2.C1 <- plot2.C1 + labs(x= element_blank()) 
plot2.C1 <- plot2.C1 + labs(y= element_blank())
plot2.C1 <- plot2.C1 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))

plot2.C2 <- ggplot(Sum.C.plt,aes(y=mean,x=bac,colour=bac))+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=0.5)+geom_point(size=3)
plot2.C2 <- plot2.C2 + geom_hline(yintercept=pred.C.plt, linetype="dashed", color = "blue3",size=1)
plot2.C2 <- plot2.C2 + theme_classic() + ggtitle("Churchill, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot2.C2 <- plot2.C2 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot2.C2 <- plot2.C2 + theme(legend.position = "none")
plot2.C2 <- plot2.C2 + labs(x= element_blank()) 
plot2.C2 <- plot2.C2 + labs(y= element_blank())
plot2.C2 <- plot2.C2 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))

#Wellspring

melt.W.abs<-melt(df.W.2020, id.vars=c("bac","plt"), measure.vars= "abs", na.rm = T)
Sum.W.abs<- ddply(melt.W.abs, c("bac","plt","variable"), summarise,
                mean = mean(value), sd = sd(value), count=n(),
                sem = sd(value)/sqrt(length(value)))

Sum.W.abs$bac <- factor(Sum.W.abs$bac, levels = c("Control","Flaviflagellibacter \ndeserti","Rhizorhabdus \nwittichii 2","Rhizorhabdus \nwittichii 1","Rhizobium \ncapsici 2","Rhizobium \ncapsici 1","Sphingomonas \npituitosa 2","Sphingomonas \npituitosa 1","Fervidobacterium \nriparium","Rhizobium \nrosettiformans","Pseudomonas \nprotogens","All 10 bacteria"))
Sum.W.abs$CI<-qnorm(0.975)*Sum.W.abs$sd/sqrt(Sum.W.abs$count)

Sum.W.plt<-subset(Sum.W.abs, Sum.W.abs$plt=="Y")
Sum.W.plt<-Sum.W.plt[2:12,]
Sum.W.plt$CI<-qnorm(0.975)*Sum.W.plt$sd/sqrt(Sum.W.plt$count)

Sum.W.bac<-subset(Sum.W.abs, Sum.W.abs$plt=="N")
Sum.W.bac$CI<-qnorm(0.975)*Sum.W.bac$sd/sqrt(Sum.W.bac$count)

pred.W.plt<-mean(Sum.W.plt$mean[1:10])
pred.W.bac<-mean(Sum.W.bac$mean[1:10])

plot2.W1 <- ggplot(Sum.W.bac,aes(y=mean,x=bac,colour=bac))+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=0.5)+geom_point(size=3)
plot2.W1 <- plot2.W1 + geom_hline(yintercept=pred.W.bac, linetype="dashed", color = "blue3",size=1)
plot2.W1 <- plot2.W1 + theme_classic() + ggtitle("Wellspring, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot2.W1 <- plot2.W1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot2.W1 <- plot2.W1 + theme(legend.position = "none")
plot2.W1 <- plot2.W1 + labs(x= element_blank()) 
plot2.W1 <- plot2.W1 + labs(y= element_blank())
plot2.W1 <- plot2.W1 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))

plot2.W2 <- ggplot(Sum.W.plt,aes(y=mean,x=bac,colour=bac))+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=0.5)+geom_point(size=3)
plot2.W2 <- plot2.W2 + geom_hline(yintercept=pred.W.plt, linetype="dashed", color = "blue3",size=1)
plot2.W2 <- plot2.W2 + theme_classic() + ggtitle("Wellspring, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot2.W2 <- plot2.W2 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot2.W2 <- plot2.W2 + theme(legend.position = "none")
plot2.W2 <- plot2.W2 + labs(x= element_blank()) 
plot2.W2 <- plot2.W2 + labs(y= element_blank())
plot2.W2 <- plot2.W2 + theme(axis.text.x = element_text(face="italic",angle = 90,vjust=0.5,hjust=1))

# Re-scale y axes
plot2.C1 <- plot2.C1 + ylim(0,11000)
plot2.W1 <- plot2.W1 + ylim(0,11000)

plot2.C2 <- plot2.C2 + ylim(0,10500)
plot2.W2 <- plot2.W2 + ylim(0,10500)

#OK, now we need to combine them both and make a reaction norm plot.

Sum.C.abs2 <- Sum.C.abs[2:23,]

plot2.C3 <- ggplot(Sum.C.abs2, aes(x = plt, y = mean, colour=bac)) +  stat_summary(aes(group = bac), fun = mean, geom = "path") +geom_point(size=3)
plot2.C3 <- plot2.C3 + theme_classic() + ggtitle("Churchill, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot2.C3 <- plot2.C3 + labs(x= element_blank())
plot2.C3 <- plot2.C3 + labs(y= element_blank())
plot2.C3 <- plot2.C3 + theme(axis.title = element_text(face="bold", size=12))
plot2.C3 <- plot2.C3 + scale_x_discrete(labels=c("N" = "No duckweed", "Y" = "Duckweed"))
plot2.C3 <- plot2.C3 + theme(axis.text.x = element_text(face="bold",size = 12, colour="black"))
plot2.C3 <- plot2.C3 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot2.C3 <- plot2.C3 + theme(legend.position = "none")

Sum.W.abs2 <- Sum.W.abs[2:23,]

plot2.W3 <- ggplot(Sum.W.abs2, aes(x = plt, y = mean, colour=bac)) +  stat_summary(aes(group = bac), fun = mean, geom = "path") +geom_point(size=3)
plot2.W3 <- plot2.W3 + theme_classic() + ggtitle("Wellspring, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot2.W3 <- plot2.W3 + labs(x= element_blank())
plot2.W3 <- plot2.W3 + labs(y= element_blank())
plot2.W3 <- plot2.W3 + theme(axis.title = element_text(face="bold", size=12))
plot2.W3 <- plot2.W3 + scale_x_discrete(labels=c("N" = "No duckweed", "Y" = "Duckweed"))
plot2.W3 <- plot2.W3 + theme(axis.text.x = element_text(face="bold",size = 12, colour="black"))
plot2.W3 <- plot2.W3 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot2.W3 <- plot2.W3 + theme(legend.position = "none")

plot2.bac1<-plot_grid(plot2.C1, plot2.C2, plot2.C3, ncol=3,nrow=1, align='h', labels = c('A','B','C'))
plot2.bac2<-plot_grid(plot2.W1, plot2.W2, plot2.W3, ncol=3,nrow=1, align='h', labels = c('D','E','F'))

plot2.1<-plot_grid(plot2.bac1,plot2.bac2, nrow = 2, ncol=1, align='v')

y.grob_abs <- textGrob(expression(bold("                       Microbial density (cells/µL)")), gp=gpar(fontsize=12), rot=90)

x.grob_bac2 <- textGrob(expression(bold("Bacteria                                                                     ")), gp=gpar(fontsize=12))

x.grob.plt <- textGrob(expression(bold("No plants                                                                             Plants                                                                                ")), gp=gpar(fontsize=12))

grid.arrange(arrangeGrob(plot2.1, left = y.grob_abs, bottom = x.grob_bac2, top=x.grob.plt))
```

![](Ch1_final_fuller_code_files/figure-html/Microbial growth models-3.png)<!-- -->

```r
#Let's arrange it the other way too.

plot2.C1b <- plot2.C1 + labs(y= "No plants") 
plot2.C1b <- plot2.C1b + theme(axis.title = element_text(size=12, face='bold'))

plot2.C2b <- plot2.C2 + labs(y= "Plants") 
plot2.C2b <- plot2.C2b + theme(axis.title = element_text(size=12, face='bold'))

plot2.bac1b<-plot_grid(plot2.C1b, plot2.W1, plot2.C3, ncol=3,nrow=1, align='h', labels = c('A','B','E'))
plot2.bac2b<-plot_grid(plot2.C2b, plot2.W2, plot2.W3, ncol=3,nrow=1, align='h', labels = c('C','D','F'))

plot2.2<-plot_grid(plot2.bac1b,plot2.bac2b, nrow = 2, ncol=1, align='v')

grid.arrange(arrangeGrob(plot2.2, left = y.grob_abs, bottom = x.grob_bac2))
```

![](Ch1_final_fuller_code_files/figure-html/Microbial growth models-4.png)<!-- -->

SECTION 5: Fitness regression


```r
# We will fit models (with full random effect structure) from which we can extract 'genotype' (i.e. bacterial inocula) means for their (1) effects on duckweed growth, and (2) microbial density when grown with plants.
# We will use scaled data so that we can center fitness proxy data around 0.

df.C.fit<-subset(df.C.2020, df.C.2020$bac!="Control" & df.C.2020$plt=="Y")
df.C.fit$Grt.scaled<-(df.C.fit$Grt-mean(df.C.fit$Grt,na.rm=T))/sd(df.C.fit$Grt,na.rm=T)
df.C.fit$pix0.scaled<-(df.C.fit$pix0-mean(df.C.fit$pix0,na.rm=T))/sd(df.C.fit$pix0,na.rm=T)
df.C.fit$abs.scaled<-(df.C.fit$abs-mean(df.C.fit$abs,na.rm=T))/sd(df.C.fit$abs,na.rm=T)

# Now we create models that we can extract strain/genotype means from
modC.grt<-lmer(Grt.scaled ~ bac + pix0.scaled + (1|edge) +(1|plate), data=df.C.fit)
em.C.grt<-as.data.frame(emmeans(modC.grt, "bac", var="Grt.scaled"))
em.C.grt$CI<-em.C.grt$SE*qnorm(0.975)

modC.abs<-lmer(abs.scaled ~ bac + (1|edge) +(1|plate), data=df.C.fit)
```

```
## boundary (singular) fit: see help('isSingular')
```

```r
em.C.abs<-as.data.frame(emmeans(modC.abs, "bac", var="abs.scaled"))
em.C.abs$CI<-em.C.abs$SE*qnorm(0.975)

em.C.grt$abs<-em.C.abs$emmean
em.C.grt$abs.CI<-em.C.abs$CI
em.C.grt$abs.SE<-em.C.abs$SE

em.C.grt$num<-1:11

# Is there a significant correlation between bacterial cell density and duckweed growth?
regC<-lm(emmean~abs,data=em.C.grt)
summary(regC)
```

```
## 
## Call:
## lm(formula = emmean ~ abs, data = em.C.grt)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.24481 -0.11200 -0.05293  0.06383  0.46403 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept) -0.07871    0.06023  -1.307  0.22367   
## abs          0.26343    0.06942   3.795  0.00425 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1988 on 9 degrees of freedom
## Multiple R-squared:  0.6154,	Adjusted R-squared:  0.5726 
## F-statistic:  14.4 on 1 and 9 DF,  p-value: 0.004252
```

```r
plot3.C1 <- ggplot(em.C.grt, aes(x = abs, y = emmean, colour = bac, fill=bac)) + geom_point(size=3,key_glyph=draw_key_blank) + geom_pointrange(aes(ymin=emmean-SE,ymax=emmean+SE),key_glyph=draw_key_blank) + geom_pointrange(aes(xmin=abs-abs.SE, xmax=abs+abs.SE),key_glyph=draw_key_blank)
plot3.C1 <- plot3.C1 + theme_classic() + ggtitle("Churchill, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot3.C1 <- plot3.C1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot3.C1 <- plot3.C1 + labs(x= element_blank())
plot3.C1 <- plot3.C1 + labs(y= "Duckweed growth")
plot3.C1 <- plot3.C1 + theme(axis.title = element_text(face="bold", size=12))
plot3.C1 <- plot3.C1 + geom_text_repel(aes(label = num), point.padding = 0.25, box.padding =0.25, fontface="bold", show.legend = FALSE)
plot3.C1 <- plot3.C1 + theme(legend.title=element_blank())
plot3.C1 <- plot3.C1 + scale_fill_discrete(labels = c(
  substitute(paste(bold("1 "), bolditalic("Flavobacterium succinicans "),bold("1"))), 
  substitute(paste(bold("2 "), bolditalic("Bosea massiliensis"))),
  substitute(paste(bold("3 "), bolditalic("Aeromonas salmonicida"))), 
  substitute(paste(bold("4 "), bolditalic("Ohtaekwangia koreensis"))), 
  substitute(paste(bold("5 "), bolditalic("Flavobacterium succinicans "),bold("2"))),
  substitute(paste(bold("6 "), bolditalic("Falsiroseomonas stagni"))), 
  substitute(paste(bold("7 "), bolditalic("Parasediminibacterium paludis"))), 
  substitute(paste(bold("8 "), bolditalic("Arcicella "), bold("sp."))), 
  substitute(paste(bold("9 "), bolditalic("Microbacterium oxydans"))), 
  substitute(paste(bold("10 "), bolditalic("Pseudomonas protogens"))),
  substitute(paste(bold("11 "), bold("All 10 bacteria")))
))
plot3.C1 <- plot3.C1 + guides(colour = FALSE)
```

```
## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
## "none")` instead.
```

```r
df.W.fit<-subset(df.W.2020, df.W.2020$bac!="Control" & df.W.2020$plt=="Y")
df.W.fit$Grt.scaled<-(df.W.fit$Grt-mean(df.W.fit$Grt,na.rm=T))/sd(df.W.fit$Grt,na.rm=T)
df.W.fit$pix0.scaled<-(df.W.fit$pix0-mean(df.W.fit$pix0,na.rm=T))/sd(df.W.fit$pix0,na.rm=T)
df.W.fit$abs.scaled<-(df.W.fit$abs-mean(df.W.fit$abs,na.rm=T))/sd(df.W.fit$abs,na.rm=T)

modW.grt<-lmer(Grt.scaled ~ bac + pix0.scaled + (1|edge) +(1|plate), data=df.W.fit)
em.W.grt<-as.data.frame(emmeans(modW.grt, "bac", var="Grt.scaled"))
em.W.grt$CI<-em.W.grt$SE*qnorm(0.975)

modW.abs<-lmer(abs.scaled ~ bac + (1|edge) +(1|plate), data=df.W.fit)
em.W.abs<-as.data.frame(emmeans(modW.abs, "bac", var="abs.scaled"))
em.W.abs$CI<-em.W.abs$SE*qnorm(0.975)

em.W.grt$abs<-em.W.abs$emmean
em.W.grt$abs.CI<-em.W.abs$CI
em.W.grt$abs.SE<-em.W.abs$SE
em.W.grt$num<-1:11

# Is there a significant correlation between bacterial cell density and duckweed growth?
regW<-lm(emmean~abs,data=em.W.grt)
summary(regW)
```

```
## 
## Call:
## lm(formula = emmean ~ abs, data = em.W.grt)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.19245 -0.10656  0.01918  0.07309  0.19960 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept) -0.05355    0.04099  -1.306   0.2238  
## abs          0.16355    0.06141   2.663   0.0259 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1347 on 9 degrees of freedom
## Multiple R-squared:  0.4408,	Adjusted R-squared:  0.3786 
## F-statistic: 7.093 on 1 and 9 DF,  p-value: 0.0259
```

```r
plot3.W1 <- ggplot(em.W.grt, aes(x = abs, y = emmean, colour = bac, fill=bac)) + geom_point(size=3,key_glyph=draw_key_blank) + geom_pointrange(aes(ymin=emmean-SE,ymax=emmean+SE),key_glyph=draw_key_blank) + geom_pointrange(aes(xmin=abs-abs.SE, xmax=abs+abs.SE),key_glyph=draw_key_blank)
plot3.W1 <- plot3.W1 + theme_classic() + ggtitle("Wellspring, 2020") + theme(plot.title = element_text(size=14, face='bold', hjust = 0.5))
plot3.W1 <- plot3.W1 + scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","blue3"))
plot3.W1 <- plot3.W1 + labs(x= element_blank())
plot3.W1 <- plot3.W1 + labs(y= element_blank())
plot3.W1 <- plot3.W1 + theme(axis.title = element_text(face="bold", size=12))
plot3.W1 <- plot3.W1 + geom_text_repel(aes(label = num), point.padding = 0.25, box.padding =0.25, fontface="bold", show.legend = FALSE)
plot3.W1 <- plot3.W1 + theme(legend.title=element_blank())
plot3.W1 <- plot3.W1 + scale_fill_discrete(labels = c(
  substitute(paste(bold("1 "), bolditalic("Sphingomonas pituitosa "), bold("1"))), 
  substitute(paste(bold("2 "), bolditalic("Flaviflagellibacter deserti"))),
  substitute(paste(bold("3 "), bolditalic("Rhizobium rosettiformans"))), 
  substitute(paste(bold("4 "), bolditalic("Rhizorhabdus wittichii "), bold("1"))), 
  substitute(paste(bold("5 "), bolditalic("Rhizobium capsici "), bold("1"))),
  substitute(paste(bold("6 "), bolditalic("Pseudomonas protogens"))), 
  substitute(paste(bold("7 "), bolditalic("Rhizorhabdus wittichii "), bold("2"))), 
  substitute(paste(bold("8 "), bolditalic("Rhizobium capsici "), bold("2"))),
  substitute(paste(bold("9 "), bolditalic("Sphingomonas pituitosa "), bold("2"))),
  substitute(paste(bold("10 "), bolditalic("Fervidobacterium riparium"))),
  substitute(paste(bold("11 "), bold("All 10 bacteria")))
))
plot3.W1 <- plot3.W1 + guides(colour = FALSE)
```

```
## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
## "none")` instead.
```

```r
# Re-scale axes?
#plot3.C1 <- plot3.C1 + ylim(-0.75,1)
#plot3.W1 <- plot3.W1 + ylim(-0.75,1)
#plot3.C1 <- plot3.C1 + xlim(-1.5,3)
#plot3.W1 <- plot3.W1 + xlim(-1.5,3)

plot3<-plot_grid(plot3.C1, plot3.W1, ncol=2,nrow=1, align='h', labels = c('A','B'))

x.grob_abs <- textGrob(expression(bold("Microbial density")), gp=gpar(fontsize=12))
grid.arrange(arrangeGrob(plot3, bottom = x.grob_abs))
```

![](Ch1_final_fuller_code_files/figure-html/Fitness regression-1.png)<!-- -->

```r
# 15 x 5?
```
