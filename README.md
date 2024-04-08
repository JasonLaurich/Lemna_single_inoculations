Single versus 10-strain inoculations of bacteria on Lemna
================
Jason Laurich and Megan Frederickson
Sys.Date()

This code generates the figures and model results for:

Laurich JR, Lash E, O’Brien AM, Pogoutse O, Frederickson ME. Community
interactions among microbes give rise to host-microbiome mutualisms in
an aquatic plant.

### SECTION 1: Upload packages

``` r
library(lme4) # mixed effects models
library(lmerTest) # get results of lmer models
library(car) #type 3 anovas
library(tidyverse) #manipulating data frames
library(cowplot) # luxury plotting
library(emmeans) # allows extraction of bacteria-specific regression estimates
library(ggrepel) # fine tune figures
library(ggtext) # fine tune figures
library(ggforce) # fine tune figures even more
library(kableExtra) # exports tables to LaTeX format
library(ape) #phylogenetics
library(phytools) #phylogenetics
library(caper) #phylogenetics
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggtree) #luxury tree plotting
#devtools::install_github('datarootsio/artyfarty')
library(artyfarty) #ggplot themes
library(magick) #graphics
```

### SECTION 2: Upload data

``` r
#Load data
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

#Make edge a factor
df.CW$edge <- as.factor(df.CW$edge)

# Create vector for plate - 38 full, 39th with 8 wells
plate<-vector()
for (i in 1:38){
  plate<-c(plate, rep(i,24))
}
plate<-c(plate,rep(39,8))
df.CW$plate <- plate

#Makes plate number a factor
df.CW$plate <- as.factor(df.CW$plate)

#Take out errors (pipetting errors, wrong bacteria)
df.CW<-df.CW[-c(453,494,798),]

# Calculate growth.
df.CW$Grt<-df.CW$pix1 - df.CW$pix0

# Log abs count
df.CW$logabs<-log10(df.CW$abs)

#Let's create separate dataframes for each population 
df.C.2020<-subset(df.CW, df.CW$pop== "Churchill")
df.W.2020<-subset(df.CW, df.CW$pop== "Wellspring")

# Assign names from sequencing
#Churchill
df.C.2020$bac_name <- ifelse(df.C.2020$bac == 1, "Flavobacterium \nsp.", 
                             ifelse(df.C.2020$bac == 2, "Pseudomonas \nprotegens 1", 
                                    ifelse(df.C.2020$bac == 3, "Devosia \nconfluentis", 
                                           ifelse(df.C.2020$bac == 4, "Arcicella \nsp.", 
                                                  ifelse(df.C.2020$bac == 5, "Microbacterium \noxydans", 
                                                         ifelse(df.C.2020$bac == 6, "Aeromonas salmonicida", 
                                                                ifelse(df.C.2020$bac == 7, "Unidentified \nFlavobacteriaceae",
                                                                       ifelse(df.C.2020$bac == 8, "Bosea \nmassiliensis",                                                                                                 ifelse(df.C.2020$bac == 9, "Unidentified \nChitinophagaceae", 
                                                                                    ifelse(df.C.2020$bac == 10, "Falsiroseomonas sp.", as.character(df.C.2020$bac)))))))))))

#Wellspring
df.W.2020$bac_name <- as.factor(ifelse(df.W.2020$bac == 1, "Allorhizobium \nsp. 2", 
                             ifelse(df.W.2020$bac == 2, "Sphingomonas \npituitosa 1", 
                                    ifelse(df.W.2020$bac == 3, "Pseudomonas \nprotegens 2", 
                                           ifelse(df.W.2020$bac == 4, "Rhizobium sp.", 
                                                  ifelse(df.W.2020$bac == 5, "Rhizorhabdus \nwittichii 1", 
                                                         ifelse(df.W.2020$bac == 6, "Allorhizobium \nsp. 1", 
                                                                ifelse(df.W.2020$bac == 7, "Sphingomonas \npituitosa 2",
                                                                       ifelse(df.W.2020$bac == 8, "Unidentified \nHyphomicrobiales",                                                                                                 ifelse(df.W.2020$bac == 9, "Rhizorhabdus \nwittichi 2", ifelse(df.W.2020$bac == 10, "Rhizobium \nrosettiformans", as.character(df.W.2020$bac))))))))))))

#Bind them back together
df2020 <- rbind(df.C.2020, df.W.2020)

#Trim white space
df2020$pop <- trimws(df2020$pop)
```

### SECTION 3: Visualize phylogenetic tree

``` r
#Import tree that was created in QIIME2
tree <- ape::read.tree("tree.nwk")

#Import taxonomy generated by Qiime2 for Churchill and Wellspring field microbiomes
taxmat <- read.delim("church_well_taxonomy.tsv")

#List tree tip labels not in taxonomy from qiime2 (i.e., the single strains)
addtips <- tree$tip.label[which(!tree$tip.label %in% taxmat$Feature.ID)]

#Add these tips to taxonomy
addtax <- data.frame("Feature.ID" = addtips, "Taxon" = NA, "Confidence" = NA)
taxmat <- rbind(taxmat, addtax)

#Split taxonomy across columns
ranks <- c("domain","phylum","class","order","family","genus","species")

#Clean up and format dataframe
taxonomy <- taxmat %>%
  mutate_at('Taxon',str_replace_all, "[a-z]__","") %>%
  separate(Taxon, sep = ';', into=ranks,remove = TRUE) %>%
  column_to_rownames(var = "Feature.ID") %>%
  as.matrix()

#Write file to disk
write.csv(taxonomy, file="merged_taxonomy.csv")
```

For the next step, we manually input the single strain taxonomy
information into csv file, and then uploaded it again.

``` r
#Read back in the edited taxonomy file
taxonomy <- read.csv("merged_taxonomy_edited.csv")

#Classify tips as field asvs or cultured isolates
taxonomy$Type <- ifelse(taxonomy$label %in% addtips, "Cultured isolate", "Field ASV")

#Label the cultured isolates by species name
taxonomy$tree_label <- ifelse(taxonomy$Type == "Cultured isolate", taxonomy$species, "")

#Read in feature table
featuretable <- read.csv("church_well_feature_table.csv")

#List tree tip labels not in feature table from qiime2 (i.e., the single strains)
addtips <- tree$tip.label[which(!tree$tip.label %in% featuretable$Feature.ID)]

#Add these tips to feature table
addtax <- data.frame("Feature.ID" = addtips, "Church.F.28" = NA, "Well.F.92" = NA)
addtax$Church.F.28 <- ifelse(substring(addtax$Feature.ID, 1, 1) == "C", 1,0)
addtax$Well.F.92 <- ifelse(substring(addtax$Feature.ID, 1, 1) == "W", 1,0)
featuretable <- rbind(featuretable, addtax)
colnames(featuretable) <- c("label", "Churchill", "Wellspring")

#Add single strains to feature table
taxonomy <- merge(taxonomy, featuretable, by="label")
taxonomy$Population <- ifelse(taxonomy$Churchill > 0 & taxonomy$Wellspring == 0, "Churchill", ifelse(taxonomy$Churchill == 0 & taxonomy$Wellspring > 0, "Wellspring", "Both"))

#Identify sisters, etc.
taxonomy$sister_clade <- sapply(taxonomy$label, getSisters, tree=tree)
taxonomy$sibs <- sapply(taxonomy$sister_clade, getDescendants, tree=tree)
taxonomy$tree_class <- ifelse(taxonomy$Type == "Cultured isolate", "", taxonomy$class)
taxonomy$tree_order <- ifelse(taxonomy$Type == "Cultured isolate", "", taxonomy$order)
taxonomy$tree_family <- ifelse(taxonomy$Type == "Cultured isolate", "", taxonomy$family)
taxonomy$tree_genus <- ifelse(taxonomy$Type == "Cultured isolate", "", taxonomy$genus)
taxonomy$tree_asvlabel <- ifelse(taxonomy$Type == "Cultured isolate", "", ifelse(taxonomy$genus != "", taxonomy$genus, taxonomy$family))
```

Use ggtree to plot the tree.

``` r
#Make type and population factors
taxonomy$Type <- as.factor(taxonomy$Type)
taxonomy$Population <- as.factor(taxonomy$Population)

p <- ggtree(tree, layout="circular") %<+% taxonomy + 
  geom_tippoint(aes(shape = Type, color = Population), size=1.5, alpha=0.8)+ 
  theme(legend.direction="horizontal", legend.position=c(0.5,0.1), legend.box = "horizontal")+
  geom_tiplab2(aes(label=trimws(tree_asvlabel)), size=1.5, offset=0.005)+
  geom_tiplab2(aes(label=trimws(tree_label)), size=1.5, offset=0.005, color = "blue")+
  guides(color = guide_legend(override.aes = list(size = 3)), shape=guide_legend(override.aes = list(size=3)))

#p <- p + geom_text(aes(label=node), size=0.8, vjust=-2) #This visualizes node numbers, if needed

pdf("Figure_2.pdf", width=8, height=9)
p
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### SECTION 3: Microbial growth models and figures

First, let’s fit a model combining all data (with bacterial strain
nested within population as a random effect) in which we test whether
the productivity of 10-strain communities differs from the productivity
of single strains. Again, this tests whether there is a
biodiversity-ecosystem function relationship, whereby the number of
strains in the microbiome (1 or 10) predicts microbial productivity.

``` r
#Group by whether 0, 1, or 10 strains were added
df2020$group <- as.factor(ifelse(df2020$bac == "All10", "10-strain community", ifelse(df2020$bac == "Control", "Control", "Single strain")))

#Adjust label names
df2020[df2020$group == "10-strain community", "bac_name"] <- "All 10 bacteria"

#Make single strain the reference
df2020 = df2020 %>% mutate(group = relevel(group, "Single strain"))

#hist(log(df2020$abs))
df2020$concat <- paste0(df2020$pop, " ", df2020$bac_name) #Concatenate pop and bac names, since there are different bacteria in each population
df2020$group <- as.factor(df2020$group)
df2020$concat <- as.factor(df2020$concat)
lmer.abs <-lmer(log(as.numeric(abs))~group*plt+(1|edge)+(1|plate)+(1|concat), data=subset(df2020, bac != "Control"))
summary(lmer.abs)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(as.numeric(abs)) ~ group * plt + (1 | edge) + (1 | plate) +  
    ##     (1 | concat)
    ##    Data: subset(df2020, bac != "Control")
    ## 
    ## REML criterion at convergence: 1294.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6264 -0.5250 -0.0247  0.6230  3.0547 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  concat   (Intercept) 0.616973 0.78548 
    ##  plate    (Intercept) 0.081136 0.28484 
    ##  edge     (Intercept) 0.008628 0.09289 
    ##  Residual             0.968392 0.98407 
    ## Number of obs: 434, groups:  concat, 22; plate, 19; edge, 2
    ## 
    ## Fixed effects:
    ##                               Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)                     5.4636     0.2112  18.3398  25.870 6.86e-16 ***
    ## group10-strain community        2.6116     0.6263  22.3735   4.170 0.000387 ***
    ## pltY                            2.1637     0.1011 404.3985  21.398  < 2e-16 ***
    ## group10-strain community:pltY  -1.6261     0.3670 397.7757  -4.431 1.21e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) gr10-c pltY  
    ## grp10-strnc -0.273              
    ## pltY        -0.238  0.083       
    ## grp10-cmm:Y  0.068 -0.228 -0.269

``` r
#Anova(lmer.abs, type=3)
ranova(lmer.abs) #Random effects
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## log(as.numeric(abs)) ~ group + plt + (1 | edge) + (1 | plate) + (1 | concat) + group:plt
    ##              npar  logLik    AIC     LRT Df Pr(>Chisq)    
    ## <none>          8 -647.15 1310.3                          
    ## (1 | edge)      7 -647.48 1309.0   0.656  1  0.4179127    
    ## (1 | plate)     7 -653.62 1321.2  12.940  1  0.0003216 ***
    ## (1 | concat)    7 -709.78 1433.5 125.255  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#plot(lmer.abs) #Looks okay
```

Generate LaTeX Table 1 for the manuscript.

``` r
tbl <- as.data.frame(coef(summary(lmer.abs)))
tbl[, 1:4] <- round(tbl[ , 1:4], 3)
colnames(tbl) <- c("Estimate", "SE", "df", "t-value", "p-value")
rownames(tbl) <- c("Intercept", "10-strain community", "Host", "10-strain community x host")
tbl$`p-value` <- ifelse(tbl$`p-value` < 0.001, "< 0.001", tbl$`p-value`)

kbl(tbl, caption = "Model results for the effect of microbial strain diversity (one versus ten strains) on microbial productivity. Intercept is a single microbial strain growing in the absence of a host. Estimates are on a log scale.", booktabs = T, format = "latex") %>% kable_styling(position = "center")
```

Next, let’s calculate the additive expectation for productivity in
10-strain communities, given the productivity of each single strain
grown separately.

``` r
#Make the 10-strain treatment the reference
df2020$bac_name <- as.factor(df2020$bac_name)
df2020 = df2020 %>% mutate(bac_name = relevel(bac_name, "All 10 bacteria"))

#Churchill 2020
#Model 
lmer.abs.C<-lmer(abs~bac_name*plt + (1|edge) + (1|plate), data=subset(df2020, bac_name != "Control" & pop=="Churchill"))
summary(lmer.abs.C)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs ~ bac_name * plt + (1 | edge) + (1 | plate)
    ##    Data: subset(df2020, bac_name != "Control" & pop == "Churchill")
    ## 
    ## REML criterion at convergence: 3729.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5520 -0.2659 -0.0436  0.1280  7.6019 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  plate    (Intercept)   96245   310.2  
    ##  edge     (Intercept)   29877   172.9  
    ##  Residual             4673921  2161.9  
    ## Number of obs: 224, groups:  plate, 19; edge, 2
    ## 
    ## Fixed effects:
    ##                                               Estimate Std. Error       df
    ## (Intercept)                                    3258.23     705.68    71.01
    ## bac_nameAeromonas salmonicida                 -2570.84     935.65   200.92
    ## bac_nameArcicella \nsp.                       -2407.99     934.97   201.82
    ## bac_nameBosea \nmassiliensis                  -3085.27    1035.73   201.64
    ## bac_nameDevosia \nconfluentis                 -2787.08     976.94   201.32
    ## bac_nameFalsiroseomonas sp.                   -3107.20     903.02   200.77
    ## bac_nameFlavobacterium \nsp.                  -3156.33     956.05   201.78
    ## bac_nameMicrobacterium \noxydans              -3197.74     934.98   201.13
    ## bac_namePseudomonas \nprotegens 1              4649.71    1000.20   199.20
    ## bac_nameUnidentified \nChitinophagaceae        -683.37    1128.84   199.50
    ## bac_nameUnidentified \nFlavobacteriaceae      -2713.76     934.61   201.05
    ## pltY                                           1828.70    1030.30   195.61
    ## bac_nameAeromonas salmonicida:pltY             -519.31    1331.92   197.86
    ## bac_nameArcicella \nsp.:pltY                     78.79    1374.16   197.73
    ## bac_nameBosea \nmassiliensis:pltY              -510.69    1435.68   200.69
    ## bac_nameDevosia \nconfluentis:pltY            -1466.68    1462.92   200.82
    ## bac_nameFalsiroseomonas sp.:pltY                399.04    1330.96   199.88
    ## bac_nameFlavobacterium \nsp.:pltY               908.95    1475.00   199.76
    ## bac_nameMicrobacterium \noxydans:pltY           322.86    1406.75   195.61
    ## bac_namePseudomonas \nprotegens 1:pltY         -411.63    1478.65   199.52
    ## bac_nameUnidentified \nChitinophagaceae:pltY    -33.79    1562.85   197.54
    ## bac_nameUnidentified \nFlavobacteriaceae:pltY   418.36    1407.16   196.32
    ##                                               t value Pr(>|t|)    
    ## (Intercept)                                     4.617 1.69e-05 ***
    ## bac_nameAeromonas salmonicida                  -2.748 0.006548 ** 
    ## bac_nameArcicella \nsp.                        -2.575 0.010725 *  
    ## bac_nameBosea \nmassiliensis                   -2.979 0.003248 ** 
    ## bac_nameDevosia \nconfluentis                  -2.853 0.004785 ** 
    ## bac_nameFalsiroseomonas sp.                    -3.441 0.000705 ***
    ## bac_nameFlavobacterium \nsp.                   -3.301 0.001138 ** 
    ## bac_nameMicrobacterium \noxydans               -3.420 0.000758 ***
    ## bac_namePseudomonas \nprotegens 1               4.649 6.07e-06 ***
    ## bac_nameUnidentified \nChitinophagaceae        -0.605 0.545618    
    ## bac_nameUnidentified \nFlavobacteriaceae       -2.904 0.004101 ** 
    ## pltY                                            1.775 0.077466 .  
    ## bac_nameAeromonas salmonicida:pltY             -0.390 0.697032    
    ## bac_nameArcicella \nsp.:pltY                    0.057 0.954337    
    ## bac_nameBosea \nmassiliensis:pltY              -0.356 0.722429    
    ## bac_nameDevosia \nconfluentis:pltY             -1.003 0.317273    
    ## bac_nameFalsiroseomonas sp.:pltY                0.300 0.764631    
    ## bac_nameFlavobacterium \nsp.:pltY               0.616 0.538438    
    ## bac_nameMicrobacterium \noxydans:pltY           0.230 0.818715    
    ## bac_namePseudomonas \nprotegens 1:pltY         -0.278 0.781010    
    ## bac_nameUnidentified \nChitinophagaceae:pltY   -0.022 0.982774    
    ## bac_nameUnidentified \nFlavobacteriaceae:pltY   0.297 0.766547    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#plot(lmer.abs.C) # Not great, but easier/more interpretable to leave this on a raw scale for estimating additive expectation, rather than log-transforming response variable

#Additive and mean predictions
C_emmeans <- as.data.frame(emmeans(lmer.abs.C, specs = ~bac_name+plt))
C_micr_mean_noplt <- mean(C_emmeans[C_emmeans$plt == "N" & C_emmeans$bac_name != "All 10 bacteria", "emmean"])
C_micr_mean_plt <- mean(C_emmeans[C_emmeans$plt == "Y" & C_emmeans$bac_name != "All 10 bacteria", "emmean"])
C_micr_sum_noplt <- sum(C_emmeans[C_emmeans$plt == "N" & C_emmeans$bac_name != "All 10 bacteria", "emmean"])
C_micr_sum_plt <- sum(C_emmeans[C_emmeans$plt == "Y" & C_emmeans$bac_name != "All 10 bacteria", "emmean"])

#Wellspring 2020

#Model 
lmer.abs.W<-lmer(abs~bac_name*plt + (1|edge) + (1|plate), data=subset(df2020, bac_name != "Control" & pop=="Wellspring"))
summary(lmer.abs.W)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs ~ bac_name * plt + (1 | edge) + (1 | plate)
    ##    Data: subset(df2020, bac_name != "Control" & pop == "Wellspring")
    ## 
    ## REML criterion at convergence: 3394.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6923 -0.4379 -0.0993  0.2676  3.3302 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  plate    (Intercept)  174840   418.1  
    ##  edge     (Intercept)  255761   505.7  
    ##  Residual             2959542  1720.3  
    ## Number of obs: 210, groups:  plate, 19; edge, 2
    ## 
    ## Fixed effects:
    ##                                               Estimate Std. Error        df
    ## (Intercept)                                   4290.429    642.371     7.826
    ## bac_nameAllorhizobium \nsp. 1                -4117.662    809.635   181.692
    ## bac_nameAllorhizobium \nsp. 2                -3307.629    761.587   179.633
    ## bac_namePseudomonas \nprotegens 2            -2760.081    846.801   183.242
    ## bac_nameRhizobium \nrosettiformans           -3339.743    765.632   181.802
    ## bac_nameRhizobium sp.                        -4150.789    788.669   185.392
    ## bac_nameRhizorhabdus \nwittichi 2            -4198.968    704.070   182.875
    ## bac_nameRhizorhabdus \nwittichii 1           -4098.765    820.393   186.515
    ## bac_nameSphingomonas \npituitosa 1           -3404.490    796.308   186.790
    ## bac_nameSphingomonas \npituitosa 2           -3720.620    762.618   182.181
    ## bac_nameUnidentified \nHyphomicrobiales      -4127.509    793.843   186.568
    ## pltY                                          2495.212    948.156   181.346
    ## bac_nameAllorhizobium \nsp. 1:pltY           -1514.928   1218.086   177.988
    ## bac_nameAllorhizobium \nsp. 2:pltY             -23.689   1267.062   185.214
    ## bac_namePseudomonas \nprotegens 2:pltY        -238.306   1281.919   185.655
    ## bac_nameRhizobium \nrosettiformans:pltY        974.289   1217.264   183.511
    ## bac_nameRhizobium sp.:pltY                    2151.108   1217.382   180.957
    ## bac_nameRhizorhabdus \nwittichi 2:pltY         799.189   1199.187   181.789
    ## bac_nameRhizorhabdus \nwittichii 1:pltY        176.455   1291.634   184.944
    ## bac_nameSphingomonas \npituitosa 1:pltY       -824.553   1253.536   186.876
    ## bac_nameSphingomonas \npituitosa 2:pltY       3789.222   1221.788   182.488
    ## bac_nameUnidentified \nHyphomicrobiales:pltY -1627.915   1299.809   182.195
    ##                                              t value Pr(>|t|)    
    ## (Intercept)                                    6.679 0.000172 ***
    ## bac_nameAllorhizobium \nsp. 1                 -5.086 9.07e-07 ***
    ## bac_nameAllorhizobium \nsp. 2                 -4.343 2.34e-05 ***
    ## bac_namePseudomonas \nprotegens 2             -3.259 0.001331 ** 
    ## bac_nameRhizobium \nrosettiformans            -4.362 2.16e-05 ***
    ## bac_nameRhizobium sp.                         -5.263 3.89e-07 ***
    ## bac_nameRhizorhabdus \nwittichi 2             -5.964 1.25e-08 ***
    ## bac_nameRhizorhabdus \nwittichii 1            -4.996 1.34e-06 ***
    ## bac_nameSphingomonas \npituitosa 1            -4.275 3.04e-05 ***
    ## bac_nameSphingomonas \npituitosa 2            -4.879 2.32e-06 ***
    ## bac_nameUnidentified \nHyphomicrobiales       -5.199 5.22e-07 ***
    ## pltY                                           2.632 0.009229 ** 
    ## bac_nameAllorhizobium \nsp. 1:pltY            -1.244 0.215247    
    ## bac_nameAllorhizobium \nsp. 2:pltY            -0.019 0.985104    
    ## bac_namePseudomonas \nprotegens 2:pltY        -0.186 0.852728    
    ## bac_nameRhizobium \nrosettiformans:pltY        0.800 0.424519    
    ## bac_nameRhizobium sp.:pltY                     1.767 0.078915 .  
    ## bac_nameRhizorhabdus \nwittichi 2:pltY         0.666 0.505973    
    ## bac_nameRhizorhabdus \nwittichii 1:pltY        0.137 0.891485    
    ## bac_nameSphingomonas \npituitosa 1:pltY       -0.658 0.511488    
    ## bac_nameSphingomonas \npituitosa 2:pltY        3.101 0.002232 ** 
    ## bac_nameUnidentified \nHyphomicrobiales:pltY  -1.252 0.212020    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Additive and mean predictions
W_emmeans <- as.data.frame(emmeans(lmer.abs.W, specs = ~bac_name+plt))
W_micr_mean_noplt <- mean(W_emmeans[W_emmeans$plt == "N" & W_emmeans$bac_name != "All 10 bacteria", "emmean"])
W_micr_mean_plt <- mean(W_emmeans[W_emmeans$plt == "Y" & W_emmeans$bac_name != "All 10 bacteria", "emmean"])
W_micr_sum_noplt <- sum(W_emmeans[W_emmeans$plt == "N" & W_emmeans$bac_name != "All 10 bacteria", "emmean"])
W_micr_sum_plt <- sum(W_emmeans[W_emmeans$plt == "Y" & W_emmeans$bac_name != "All 10 bacteria", "emmean"])
```

Next, let’s make a figure of the microbial data.

``` r
#Set figure defaults
pt_size = 3
ylimits = c(0, 50000)
zero_line <- "red3"
sum_line <- "green3"
mean_line <- "blue3"
dot_colors <- c("blue", "red", "black")

# Generate Figure 3
Sum_2020_C_micr <- subset(df2020, pop == "Churchill") %>% group_by(bac_name, plt) %>% summarize(n=n(), mean=mean(abs, na.rm=TRUE),
sd = sd(abs, na.rm=TRUE), sem = sd/sqrt(n))

#Add colors
Sum_2020_C_micr$color <- ifelse(Sum_2020_C_micr$bac_name == "All10", "10-strain community", ifelse(Sum_2020_C_micr$bac_name == "Control", "Control", "Single strain"))

#Churchill WITHOUT PLANTS
plot_C2020_micr_noplt <- ggplot(subset(Sum_2020_C_micr, bac_name != "Control" & plt == "N"), aes(y=mean,x=reorder(bac_name, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac_name, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
  theme_cowplot()+
  xlab("Bacterial treatment")+
  ggtitle("Churchill, no plants")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(legend.title = element_blank())+
  scale_y_continuous(limits=ylimits)+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.1, 0.9))+
  geom_hline(yintercept=C_micr_mean_noplt, linetype="dashed", color = mean_line,linewidth=1)+
  geom_hline(yintercept=C_micr_sum_noplt, linetype="dashed", color = sum_line,linewidth=1)+
  ylab(expression("Microbial density (cells/µL)"))

plot_C2020_micr_noplt_updated <- plot_C2020_micr_noplt+
   scale_x_discrete(labels=c("*Microbacterium oxydans*", "*Bosea massiliensis*",  "*Flavobacterium* sp.", "*Falsiroseomonas* sp.", "*Devosia confluentis*", "Unidentified *Flavobacteriaceae*", "*Aeromonas salmonicida*", "*Arcicella* sp.", "Unidentified *Chitinophagaceae*", "All 10 bacteria", "*Pseudomonas protogens* 1"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_C2020_micr_noplt_updated_v2 <-plot_C2020_micr_noplt_updated+facet_zoom(ylim = c(0, 10000))

#Churchill WITH PLANTS
plot_C2020_micr_plt <- ggplot(subset(Sum_2020_C_micr, bac_name != "Control" & plt == "Y"), aes(y=mean,x=reorder(bac_name, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac_name, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
  theme_cowplot()+
  xlab("Bacterial treatment")+
  ggtitle("Churchill, with plants")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(legend.title = element_blank())+
  scale_y_continuous(limits=ylimits)+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.1, 0.9))+
  geom_hline(yintercept=C_micr_mean_plt, linetype="dashed", color = mean_line,linewidth=1)+
  geom_hline(yintercept=C_micr_sum_plt, linetype="dashed", color = sum_line,linewidth=1)+
  ylab(expression("Microbial density (cells/µL)"))

plot_C2020_micr_plt_updated <- plot_C2020_micr_plt+
   scale_x_discrete(labels=c("*Devosia confluentis*", "*Bosea massiliensis*", "*Aeromonas salmonicida*", "*Microbacterium oxydans*",  "*Arcicella* sp.", "*Falsiroseomonas* sp.", "*Flavobacterium* sp.",  "Unidentified *Flavobacteriaceae*", "Unidentified *Chitinophagaceae*", "All 10 bacteria", "*Pseudomonas protogens* 1"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_C2020_micr_plt_updated_v2 <-plot_C2020_micr_plt_updated+facet_zoom(ylim = c(0, 10000))

# Wellspring
Sum_2020_W_micr <- subset(df2020, pop=="Wellspring") %>% group_by(bac_name, plt) %>% summarize(n=n(), mean=mean(abs, na.rm=TRUE),
sd = sd(abs, na.rm=TRUE), sem = sd/sqrt(n))

#Add colors
Sum_2020_W_micr$color <- ifelse(Sum_2020_W_micr$bac_name == "All10", "10-strain community", ifelse(Sum_2020_W_micr$bac_name == "Control", "Control", "Single strain"))

#Wellspring WITHOUT PLANTS
plot_W2020_micr_noplt <- ggplot(subset(Sum_2020_W_micr, bac_name != "Control" & plt == "N"), aes(y=mean,x=reorder(bac_name, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac_name, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
  theme_cowplot()+
  xlab("Bacterial treatment")+
  ggtitle("Wellspring, no plants")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(legend.title = element_blank())+
  scale_y_continuous(limits=ylimits)+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.1, 0.9))+
  geom_hline(yintercept=W_micr_mean_noplt, linetype="dashed", color = mean_line,linewidth=1)+
  geom_hline(yintercept=W_micr_sum_noplt, linetype="dashed", color = sum_line,linewidth=1)+
  ylab(expression("Microbial density (cells/µL)"))

plot_W2020_micr_noplt_updated <- plot_W2020_micr_noplt+
   scale_x_discrete(labels=c("Unidentified *Hyphomicrobiales*", "*Rhizorhabdus wittichii* 1", "*Rhizorhabdus wittichii* 2", "*Rhizobium* sp.", "*Allorhizobium* sp. 1",  "*Sphingomonas pituitosa* 2", "*Sphingomonas pituitosa* 1", "*Rhizobium rosettiformans*", "*Allorhizobium* sp. 2", "*Pseudomonas protogens* 2", "All 10 bacteria"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_W2020_micr_noplt_updated_v2 <-plot_W2020_micr_noplt_updated+facet_zoom(ylim = c(0, 8000))

#Wellspring WITH PLANTS
plot_W2020_micr_plt <- ggplot(subset(Sum_2020_W_micr, bac_name != "Control" & plt == "Y"), aes(y=mean,x=reorder(bac_name, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac_name, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
  theme_cowplot()+
  xlab("Bacterial treatment")+
  ggtitle("Wellspring, with plants")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(legend.title = element_blank())+
  scale_y_continuous(limits=ylimits)+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.1, 0.9))+
  geom_hline(yintercept=W_micr_mean_plt, linetype="dashed", color = mean_line,linewidth=1)+
  geom_hline(yintercept=W_micr_sum_plt, linetype="dashed", color = sum_line,linewidth=1)+
  ylab(expression("Microbial density (cells/µL)"))

plot_W2020_micr_plt_updated <- plot_W2020_micr_plt+
   scale_x_discrete(labels=c("*Allorhizobium* sp. 1", "Unidentified *Hyphomicrobiales*", "*Sphingomonas pituitosa* 1", "*Rhizorhabdus wittichii* 1", "*Allorhizobium* sp. 2", "*Rhizorhabdus wittichii* 2", "*Pseudomonas protogens* 2", "*Rhizobium rosettiformans*", "*Rhizobium* sp.", "All 10 bacteria", "*Sphingomonas pituitosa* 2"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_W2020_micr_plt_updated_v2 <-plot_W2020_micr_plt_updated+facet_zoom(ylim = c(0, 8000))

final_microbes <- plot_grid(plot_C2020_micr_noplt_updated_v2, plot_C2020_micr_plt_updated_v2,plot_W2020_micr_noplt_updated_v2, plot_W2020_micr_plt_updated_v2, nrow=2, labels="AUTO")
save_plot("Figure_3.pdf", final_microbes, base_height=12, base_width=18)
final_microbes
```

![](README_files/figure-gfm/Microbial%20growth%20figures-1.png)<!-- -->

Next, we want to plot the strain and communities means without a host
(x-axis) versus with a host (y-axis).

``` r
Sum_2020_C_micr_wide <- Sum_2020_C_micr %>% pivot_wider(names_from = plt, values_from = n:sem)

lm_micr_C <- lm(log(mean_Y)~log(mean_N), data=Sum_2020_C_micr_wide)
summary(lm_micr_C)
```

    ## 
    ## Call:
    ## lm(formula = log(mean_Y) ~ log(mean_N), data = Sum_2020_C_micr_wide)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.16532 -0.07929  0.05828  0.26874  0.48781 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   6.0359     0.6843   8.821 1.01e-05 ***
    ## log(mean_N)   0.2939     0.1045   2.814   0.0203 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4887 on 9 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.468,  Adjusted R-squared:  0.4089 
    ## F-statistic: 7.916 on 1 and 9 DF,  p-value: 0.02026

``` r
C_alt <-ggplot(subset(Sum_2020_C_micr_wide, bac_name != "Control"))+geom_point(aes(x=mean_N, y=mean_Y, color=color))+geom_abline(aes(intercept=0, slope=1), linetype="dotted")+theme_cowplot()+xlab(expression("Microbial density without hosts (cells/µL)"))+ylab(expression("Microbial density with hosts (cells/µL)"))+
  geom_errorbar(aes(x=mean_N, ymax=mean_Y+sem_Y, ymin=mean_Y-sem_Y, color=color))+
    geom_errorbarh(aes(xmax=mean_N+sem_N, xmin=mean_N-sem_N, y=mean_Y, color=color))+
 theme(legend.title = element_blank())+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.55, 0.15))+
    ggtitle("Churchill")+
    geom_text_repel(data=subset(Sum_2020_C_micr_wide, bac_name != "Control" & bac_name !="All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac_name))),fontface = "italic", size=2.5)+
    geom_text_repel(data=subset(Sum_2020_C_micr_wide, bac_name == "All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac_name))), size=2.5)+
    scale_y_log10(limits=c(10,10000))+  
  scale_x_log10(limits=c(10,10000))+
  geom_smooth(aes(x=mean_N, y=mean_Y), method="lm", se=FALSE)
C_alt
```

![](README_files/figure-gfm/Host%20versus%20no%20host-1.png)<!-- -->

``` r
Sum_2020_W_micr_wide <- Sum_2020_W_micr %>% pivot_wider(names_from = plt, values_from = n:sem)

lm_micr_W <- lm(log(mean_Y)~log(mean_N), data=Sum_2020_W_micr_wide)
summary(lm_micr_W)
```

    ## 
    ## Call:
    ## lm(formula = log(mean_Y) ~ log(mean_N), data = Sum_2020_W_micr_wide)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8852 -0.3030  0.0069  0.2509  0.7669 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   6.1026     0.8010   7.619 3.26e-05 ***
    ## log(mean_N)   0.3179     0.1295   2.455   0.0364 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5106 on 9 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.4012, Adjusted R-squared:  0.3346 
    ## F-statistic: 6.029 on 1 and 9 DF,  p-value: 0.03643

``` r
W_alt <-ggplot(subset(Sum_2020_W_micr_wide, bac_name != "Control"))+geom_point(aes(x=mean_N, y=mean_Y, color=color))+geom_abline(aes(intercept=0, slope=1), linetype="dotted")+theme_cowplot()+xlab(expression("Microbial density without hosts (cells/µL)"))+ylab(expression("Microbial density with hosts (cells/µL)"))+
 theme(legend.title = element_blank())+
    geom_errorbar(aes(x=mean_N, ymax=mean_Y+sem_Y, ymin=mean_Y-sem_Y, color=color))+
    geom_errorbarh(aes(xmax=mean_N+sem_N, xmin=mean_N-sem_N, y=mean_Y, color=color))+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.55, 0.15))+
    ggtitle("Wellspring")+
  geom_text_repel(data=subset(Sum_2020_W_micr_wide, bac_name != "Control" & bac_name !="All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac_name))),fontface = "italic", size=2.5)+
    geom_text_repel(data=subset(Sum_2020_W_micr_wide, bac_name == "All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac_name))), size=2.5)+
    scale_y_log10(limits=c(10,10000))+  
  scale_x_log10(limits=c(10,10000))+
  geom_smooth(aes(x=mean_N, y=mean_Y), method="lm", se=FALSE)
W_alt
```

![](README_files/figure-gfm/Host%20versus%20no%20host-2.png)<!-- -->

``` r
alt <- plot_grid(C_alt, W_alt, labels="AUTO")
alt
```

![](README_files/figure-gfm/Host%20versus%20no%20host-3.png)<!-- -->

``` r
save_plot("Figure_4.pdf", alt, base_height=6, base_width=12)
```

### SECTION 4: Plant growth models and figures

The first set of models tests whether there is a biodiversity-ecosystem
function relationship, whereby the number of strains in the microbiome
(0, 1, or 10) predicts plant growth.

``` r
#Make the control (uninoculated) treatment the reference
df2020 = df2020 %>% mutate(group = relevel(group, "Single strain"))
df2020 = df2020 %>% mutate(group = relevel(group, "Control"))

#Model both populations together, with bacterial treatment nested in population as a random effect
#Frond number
#hist(df2020$frd1)
lm_both_frd <- lmer(frd1~frd0+group+(1|edge)+(1|plate)+(1|pop/bac_name), data=df2020)
summary(lm_both_frd)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: frd1 ~ frd0 + group + (1 | edge) + (1 | plate) + (1 | pop/bac_name)
    ##    Data: df2020
    ## 
    ## REML criterion at convergence: 2377.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0495 -0.6193  0.0055  0.6886  3.4129 
    ## 
    ## Random effects:
    ##  Groups       Name        Variance Std.Dev.
    ##  plate        (Intercept) 0.4465   0.6682  
    ##  bac_name:pop (Intercept) 0.0377   0.1942  
    ##  pop          (Intercept) 0.4879   0.6985  
    ##  edge         (Intercept) 0.2628   0.5126  
    ##  Residual                 8.3124   2.8831  
    ## Number of obs: 475, groups:  plate, 39; bac_name:pop, 24; pop, 2; edge, 2
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)                4.8002     0.9011   6.9386   5.327  0.00112 ** 
    ## frd0                       2.0736     0.2395 464.1436   8.657  < 2e-16 ***
    ## groupSingle strain         0.4275     0.5077  19.6914   0.842  0.40983    
    ## group10-strain community   1.7780     0.6919  20.5258   2.570  0.01806 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) frd0   grpSns
    ## frd0        -0.482              
    ## grpSnglstrn -0.512 -0.001       
    ## grp10-strnc -0.391  0.031  0.673

``` r
Anova(lm_both_frd, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: frd1
    ##               Chisq Df Pr(>Chisq)    
    ## (Intercept) 28.3776  1  9.981e-08 ***
    ## frd0        74.9399  1  < 2.2e-16 ***
    ## group        8.0415  2    0.01794 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(lm_both_frd, list(pairwise ~ group), adjust="tukey")
```

    ## $`emmeans of group`
    ##  group               emmean    SE   df lower.CL upper.CL
    ##  Control               8.54 0.791 4.09     6.36     10.7
    ##  Single strain         8.97 0.643 1.83     5.94     12.0
    ##  10-strain community  10.32 0.800 4.21     8.14     12.5
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of group`
    ##  1                                     estimate    SE   df t.ratio p.value
    ##  Control - Single strain                 -0.428 0.509 19.7  -0.840  0.6833
    ##  Control - (10-strain community)         -1.778 0.695 20.5  -2.558  0.0467
    ##  Single strain - (10-strain community)   -1.350 0.516 20.6  -2.618  0.0412
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
#plot(lm_both_frd) #Looks ok

#Pixel area
#hist(df2020$pix1)
lm_both_pix <- lmer(pix1~scale(pix0)+group+(1|edge)+(1|plate)+(1|pop/bac_name), data=df2020)
summary(lm_both_pix)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pix1 ~ scale(pix0) + group + (1 | edge) + (1 | plate) + (1 |  
    ##     pop/bac_name)
    ##    Data: df2020
    ## 
    ## REML criterion at convergence: 9985.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.76518 -0.54633  0.05158  0.62057  2.87974 
    ## 
    ## Random effects:
    ##  Groups       Name        Variance Std.Dev.
    ##  plate        (Intercept) 10596583 3255    
    ##  bac_name:pop (Intercept)   379421  616    
    ##  pop          (Intercept) 31211567 5587    
    ##  edge         (Intercept)  8116286 2849    
    ##  Residual                 82393701 9077    
    ## Number of obs: 475, groups:  plate, 39; bac_name:pop, 24; pop, 2; edge, 2
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)              22258.34    4725.52     1.86   4.710   0.0485 *  
    ## scale(pix0)               6756.27     433.57   456.84  15.583   <2e-16 ***
    ## groupSingle strain        1362.99    1614.19    20.01   0.844   0.4084    
    ## group10-strain community  5752.16    2203.46    20.99   2.611   0.0163 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) scl(0) grpSns
    ## scale(pix0) -0.003              
    ## grpSnglstrn -0.311  0.007       
    ## grp10-strnc -0.229  0.030  0.676

``` r
Anova(lm_both_pix, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: pix1
    ##                Chisq Df Pr(>Chisq)    
    ## (Intercept)  22.1863  1  2.474e-06 ***
    ## scale(pix0) 242.8242  1  < 2.2e-16 ***
    ## group         8.3724  2     0.0152 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(lm_both_pix)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## pix1 ~ scale(pix0) + group + (1 | edge) + (1 | plate) + (1 | bac_name:pop) + (1 | pop)
    ##                    npar  logLik   AIC     LRT Df Pr(>Chisq)    
    ## <none>                9 -4992.8 10004                          
    ## (1 | edge)            8 -5000.8 10018 16.0662  1  6.117e-05 ***
    ## (1 | plate)           8 -5002.3 10021 19.0147  1  1.297e-05 ***
    ## (1 | bac_name:pop)    8 -4992.8 10002  0.0683  1     0.7938    
    ## (1 | pop)             8 -5006.9 10030 28.3131  1  1.032e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(lm_both_pix, list(pairwise ~ group), adjust="tukey")
```

    ## $`emmeans of group`
    ##  group               emmean   SE   df lower.CL upper.CL
    ##  Control              22242 4727 1.86      357    44126
    ##  Single strain        23605 4494 1.52    -2909    50118
    ##  10-strain community  27994 4738 1.88     6248    49739
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of group`
    ##  1                                     estimate   SE   df t.ratio p.value
    ##  Control - Single strain                  -1363 1617 19.7  -0.843  0.6815
    ##  Control - (10-strain community)          -5752 2209 20.7  -2.604  0.0424
    ##  Single strain - (10-strain community)    -4389 1633 20.4  -2.688  0.0357
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
#plot(lm_both_pix) #Looks ok
```

The second set of models fits each population separately and each
bacterium as a fixed effect to calculate additive expectation for
10-strain community from single-strain inoculations.

``` r
#Make the control treatment the reference
df2020 = df2020 %>% mutate(bac_name = relevel(bac_name, "Control"))

#Churchill
C_pixels <- lmer(pix1~pix0+bac_name+(1|edge)+(1|plate), data=subset(df2020, pop == "Churchill"))
summary(C_pixels)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pix1 ~ pix0 + bac_name + (1 | edge) + (1 | plate)
    ##    Data: subset(df2020, pop == "Churchill")
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
    ##                                           Estimate Std. Error        df t value
    ## (Intercept)                              3878.8193  3030.8194    4.8295   1.280
    ## pix0                                        4.9076     0.5278  220.0131   9.298
    ## bac_nameAll 10 bacteria                  7874.6934  2406.2248  214.3893   3.273
    ## bac_nameAeromonas salmonicida            -182.9379  2352.5254  217.5792  -0.078
    ## bac_nameArcicella \nsp.                  2258.5016  2364.4514  212.8589   0.955
    ## bac_nameBosea \nmassiliensis             1614.9724  2355.4880  216.7209   0.686
    ## bac_nameDevosia \nconfluentis            -558.6037  2337.8576  213.4642  -0.239
    ## bac_nameFalsiroseomonas sp.               707.6477  2310.0205  207.2307   0.306
    ## bac_nameFlavobacterium \nsp.             -313.0947  2380.2982  211.3749  -0.132
    ## bac_nameMicrobacterium \noxydans         1791.5725  2333.9711  214.0507   0.768
    ## bac_namePseudomonas \nprotegens 1        6424.0758  2344.7823  213.4976   2.740
    ## bac_nameUnidentified \nChitinophagaceae  2265.6999  2314.7662  206.6954   0.979
    ## bac_nameUnidentified \nFlavobacteriaceae 1814.4222  2357.8891  216.2553   0.770
    ##                                          Pr(>|t|)    
    ## (Intercept)                               0.25865    
    ## pix0                                      < 2e-16 ***
    ## bac_nameAll 10 bacteria                   0.00124 ** 
    ## bac_nameAeromonas salmonicida             0.93809    
    ## bac_nameArcicella \nsp.                   0.34056    
    ## bac_nameBosea \nmassiliensis              0.49368    
    ## bac_nameDevosia \nconfluentis             0.81138    
    ## bac_nameFalsiroseomonas sp.               0.75965    
    ## bac_nameFlavobacterium \nsp.              0.89548    
    ## bac_nameMicrobacterium \noxydans          0.44357    
    ## bac_namePseudomonas \nprotegens 1         0.00667 ** 
    ## bac_nameUnidentified \nChitinophagaceae   0.32882    
    ## bac_nameUnidentified \nFlavobacteriaceae  0.44243    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#plot(C_pixels) #Looks ok

#Calculate additive and mean predictions 
C_pix_emmeans <- as.data.frame(emmeans(C_pixels, specs=~bac_name))
pred.C.2020.avg.emmeans <- mean(as.numeric(C_pix_emmeans[3:12, "emmean"])) 
pred.C.2020.sum.emmeans <- sum(as.numeric(C_pix_emmeans[3:12, "emmean"]-as.numeric(C_pix_emmeans[1, "emmean"])))+as.numeric(C_pix_emmeans[1, "emmean"]) #Add the difference between the mean for each single strain and the control, and then add the control

#Wellspring
W_pixels <- lmer(pix1~pix0+bac_name+(1|edge)+(1|plate), data=subset(df2020, pop == "Wellspring"))
summary(W_pixels)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pix1 ~ pix0 + bac_name + (1 | edge) + (1 | plate)
    ##    Data: subset(df2020, pop == "Wellspring")
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
    ##                                           Estimate Std. Error         df
    ## (Intercept)                              5044.0552  3623.6227    10.5412
    ## pix0                                        7.2913     0.6512   223.2800
    ## bac_nameAll 10 bacteria                  4497.5803  3521.5576   222.8469
    ## bac_nameAllorhizobium \nsp. 1            1307.6764  3518.4586   223.0462
    ## bac_nameAllorhizobium \nsp. 2             659.1365  3497.7890   219.4910
    ## bac_namePseudomonas \nprotegens 2        4856.7650  3527.2609   223.3201
    ## bac_nameRhizobium \nrosettiformans       1449.3299  3580.4124   224.0491
    ## bac_nameRhizobium sp.                    4053.7219  3500.9568   219.8402
    ## bac_nameRhizorhabdus \nwittichi 2        -581.6311  3490.8712   217.3165
    ## bac_nameRhizorhabdus \nwittichii 1       3747.2298  3502.0644   220.2646
    ## bac_nameSphingomonas \npituitosa 1       1484.5091  3560.5816   220.9672
    ## bac_nameSphingomonas \npituitosa 2       5636.4924  3502.5996   220.2029
    ## bac_nameUnidentified \nHyphomicrobiales -1077.4527  3553.0940   221.5627
    ##                                         t value Pr(>|t|)    
    ## (Intercept)                               1.392    0.193    
    ## pix0                                     11.196   <2e-16 ***
    ## bac_nameAll 10 bacteria                   1.277    0.203    
    ## bac_nameAllorhizobium \nsp. 1             0.372    0.710    
    ## bac_nameAllorhizobium \nsp. 2             0.188    0.851    
    ## bac_namePseudomonas \nprotegens 2         1.377    0.170    
    ## bac_nameRhizobium \nrosettiformans        0.405    0.686    
    ## bac_nameRhizobium sp.                     1.158    0.248    
    ## bac_nameRhizorhabdus \nwittichi 2        -0.167    0.868    
    ## bac_nameRhizorhabdus \nwittichii 1        1.070    0.286    
    ## bac_nameSphingomonas \npituitosa 1        0.417    0.677    
    ## bac_nameSphingomonas \npituitosa 2        1.609    0.109    
    ## bac_nameUnidentified \nHyphomicrobiales  -0.303    0.762    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling

``` r
#plot(W_pixels) #Looks ok

#Calculate additive and mean predictions 
#emmeans
W_pix_emmeans <- as.data.frame(emmeans(W_pixels, specs=~bac_name))
pred.W.2020.avg.emmeans <- mean(as.numeric(W_pix_emmeans[3:12, "emmean"]))
pred.W.2020.sum.emmeans <- sum(as.numeric(W_pix_emmeans[3:12, "emmean"]-as.numeric(W_pix_emmeans[1, "emmean"])))+as.numeric(W_pix_emmeans[1, "emmean"]) 
```

Next, let’s make a figure of the plant growth data.

``` r
pt_size = 3
ylimits = c(12000, 48000)
zero_line <- "red3"
sum_line <- "green3"
mean_line <- "blue3"
dot_colors <- c("red","black", "blue")

Sum_2020 <- df2020 %>% group_by(pop, bac_name, group) %>% summarize(n=n(), mean=mean(pix1, na.rm=TRUE), sd=sd(pix1, na.rm=TRUE), sem=sd/sqrt(n))

plot_C2020 <- ggplot(subset(Sum_2020, pop == "Churchill"), aes(y=mean,x=reorder(bac_name, mean), color=group))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac_name,mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
  theme_cowplot()+
  geom_hline(yintercept=pred.C.2020.avg.emmeans, linetype="dashed", color = mean_line,linewidth=1)+
  geom_hline(yintercept=pred.C.2020.sum.emmeans, linetype="dashed", color = sum_line,linewidth=1)+
  xlab("Bacterial treatment")+
  ggtitle("Churchill")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
   theme(legend.title = element_blank())+
  scale_y_continuous(limits=ylimits)+
  scale_color_manual(values=dot_colors)+
  theme(legend.position = c(0.1, 0.9))+
  ylab("Duckweed area (final pixels)")

plot_C2020_updated <- plot_C2020+scale_x_discrete(labels=c("*Flavobacterium* sp.", "*Bosea massiliensis*", "Control", "Unidentified *Flavobacteriaceae*", "*Aeromonas salmonicida*", "*Devosia confluentis*", "*Falsiroseomonas* sp.", "Unidentified *Chitinophagaceae*", "*Arcicella* sp.", "*Microbacterium oxydans*", "All 10 bacteria",  "*Pseudomonas protegens*"))+
  theme(axis.text.x = ggtext::element_markdown())
plot_C2020_updated
```

![](README_files/figure-gfm/Duckweed%20growth%20figures-1.png)<!-- -->

``` r
#Wellspring 2020
plot_W2020 <- ggplot(subset(Sum_2020, pop == "Wellspring"), aes(y=mean,x=reorder(bac_name, mean), color=group))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac_name,mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
  theme_cowplot()+
  geom_hline(yintercept=pred.W.2020.avg.emmeans, linetype="dashed", color = mean_line,linewidth=1)+
  geom_hline(yintercept=pred.W.2020.sum.emmeans, linetype="dashed", color = sum_line, linewidth=1)+
  xlab("Bacterial treatment")+
  ggtitle("Wellspring")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(legend.title = element_blank())+
  scale_y_continuous(limits=ylimits)+
  scale_color_manual(values=dot_colors)+
  theme(legend.position="none")+
   ylab("Duckweed area (final pixels)")

plot_W2020_updated <- plot_W2020+scale_x_discrete(labels=c("*Sphingomonas pituitosa* 1", "Control",  "Unidentified *Hyphomicrobiales*", "*Allorhizobium* sp. 2", "*Allorhizobium* sp. 1", "*Rhizorhabdus wittichii* 2", "*Pseudomonas protegens* 2", "*Rhizorhabdus wittichii* 1", "*Rhizobium* sp.", "*Sphingomonas pituitosa* 2", "All 10 bacteria", "*Rhizobium rosettiformans*"))+
  theme(axis.text.x = ggtext::element_markdown())
plot_W2020_updated
```

![](README_files/figure-gfm/Duckweed%20growth%20figures-2.png)<!-- -->

``` r
plot_2020 <- plot_grid(plot_C2020_updated, plot_W2020_updated, labels="AUTO", align='h')
plot_2020
```

![](README_files/figure-gfm/Duckweed%20growth%20figures-3.png)<!-- -->

``` r
save_plot("Figure_5.pdf", plot_2020, base_height=8, base_width=12)
```

### SECTION 5: Fitness regression

Fit models (with full random effect structure) from which we can extract
‘genotype’ (i.e. bacterial inocula) means for their (1) effects on
duckweed growth, and (2) microbial density when grown with plants. Then,
we will use scaled data so that we can center fitness proxy data around
0.

``` r
#Churchill
df.C.2020 <- subset(df2020, bac_name != "Control" & plt == "Y" & pop == "Churchill")
df.C.2020$pix1_scaled <- scale(df.C.2020$pix1, center = TRUE, scale = TRUE)
df.C.2020$pix0_scaled <- scale(df.C.2020$pix0, center = TRUE, scale = TRUE)
df.C.2020$abs_scaled <- scale(df.C.2020$abs, center=TRUE, scale=TRUE)

# Create models to extract strain/genotype means
modC <-lmer(pix1_scaled ~ bac_name + pix0_scaled + (1|edge) +(1|plate), data=df.C.2020)
emC <- as.data.frame(emmeans(modC, "bac_name", var="pix1_scaled"))
emC$CI<-emC$SE*qnorm(0.975)

modC.abs<-lmer(abs_scaled ~ bac_name + (1|edge) +(1|plate), data=df.C.2020)
em.C.abs<-as.data.frame(emmeans(modC.abs, "bac_name", var="abs"))

emC$color <- ifelse(emC$bac_name == "All 10 bacteria", "10-strain community", "Single strain")

emC$abs<-em.C.abs$emmean
emC$abs.SE<-em.C.abs$SE

# Is there a significant correlation between bacterial cell density and duckweed growth?
regC<-lm(emmean~abs,data=emC)
summary(regC)
```

    ## 
    ## Call:
    ## lm(formula = emmean ~ abs, data = emC)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.23390 -0.10701 -0.05057  0.06099  0.44336 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept) -0.07520    0.05755  -1.307  0.22367   
    ## abs          0.25170    0.06633   3.795  0.00425 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.19 on 9 degrees of freedom
    ## Multiple R-squared:  0.6154, Adjusted R-squared:  0.5726 
    ## F-statistic:  14.4 on 1 and 9 DF,  p-value: 0.004252

``` r
plot_C_fitness <- ggplot(emC)+
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           alpha = .2)+
  geom_vline(xintercept=0, linetype="dotted")+
  geom_hline(yintercept=0, linetype="dotted")+
  geom_point(aes(x=abs, y=emmean, color=color))+theme_cowplot()+xlab(expression("Microbial density (cells/µL, centered and scaled)"))+ylab(expression("Duckweed area (pixels, centered and scaled)"))+
#geom_errorbar(aes(x=abs, ymax=emmean+SE, ymin=emmean-SE, color=color))+
#geom_errorbarh(aes(xmax=abs+abs.SE, xmin=abs-abs.SE, y=emmean, color=color))+
 theme(legend.title = element_blank())+
  scale_color_manual(values=c(dot_colors[3], dot_colors[2]))+
  theme(legend.position = c(0.55, 0.15))+
  geom_text_repel(data=subset(emC, bac_name != "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac_name))),fontface = "italic", size=2.5)+
      geom_text_repel(data=subset(emC, bac_name == "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac_name))), size=2.5)+
  geom_smooth(aes(x=abs, y=emmean), method="lm", se=FALSE)+
    scale_y_continuous(limits=c(-0.6,0.6))+
  scale_x_continuous(limits=c(-2,2.5))+
      ggtitle("Churchill")
plot_C_fitness
```

![](README_files/figure-gfm/Fitness%20regression-1.png)<!-- -->

``` r
#Wellspring
df.W.2020 <- subset(df2020, bac_name != "Control" & plt == "Y" & pop == "Wellspring")
df.W.2020$pix1_scaled <- scale(df.W.2020$pix1, center = TRUE, scale = TRUE)
df.W.2020$pix0_scaled <- scale(df.W.2020$pix0, center = TRUE, scale = TRUE)
df.W.2020$abs_scaled <- scale(df.W.2020$abs, center=TRUE, scale=TRUE)

# Create models to extract strain/genotype means
modW <-lmer(pix1_scaled ~ bac_name + pix0_scaled + (1|edge) +(1|plate), data=df.W.2020)
emW <- as.data.frame(emmeans(modW, "bac_name", var="pix1_scaled"))
emW$CI<-emW$SE*qnorm(0.975)

modW.abs<-lmer(abs_scaled ~ bac_name + (1|edge) +(1|plate), data=df.W.2020)
em.W.abs<-as.data.frame(emmeans(modW.abs, "bac_name", var="abs"))

emW$color <- ifelse(emW$bac_name == "All 10 bacteria", "10-strain community", "Single strain")

emW$abs<-em.W.abs$emmean
emW$abs.SE<-em.W.abs$SE

# Is there a significant correlation between bacterial cell density and duckweed growth?
regW<-lm(emmean~abs,data=emW)
summary(regW)
```

    ## 
    ## Call:
    ## lm(formula = emmean ~ abs, data = emW)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.18322 -0.10145  0.01826  0.06958  0.19002 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.05098    0.03903  -1.306   0.2238  
    ## abs          0.15571    0.05846   2.663   0.0259 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1283 on 9 degrees of freedom
    ## Multiple R-squared:  0.4408, Adjusted R-squared:  0.3786 
    ## F-statistic: 7.094 on 1 and 9 DF,  p-value: 0.0259

``` r
plot_W_fitness <- ggplot(emW)+annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           alpha = .2)+
  geom_vline(xintercept=0, linetype="dotted")+
  geom_hline(yintercept=0, linetype="dotted")+
geom_point(aes(x=abs, y=emmean, color=color))+theme_cowplot()+xlab(expression("Microbial density (cells/µL, centered and scaled)"))+ylab(expression("Duckweed area (pixels, centered and scaled)"))+
#geom_errorbar(aes(x=abs, ymax=emmean+SE, ymin=emmean-SE, color=color))+
#geom_errorbarh(aes(xmax=abs+abs.SE, xmin=abs-abs.SE, y=emmean, color=color))+
 theme(legend.title = element_blank())+
  scale_color_manual(values=c(dot_colors[3], dot_colors[2]))+
  theme(legend.position = c(0.55, 0.15))+
  geom_text_repel(data=subset(emW, bac_name != "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac_name))),fontface = "italic", size=2.5)+
      geom_text_repel(data=subset(emW, bac_name == "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac_name))), size=2.5)+
  geom_smooth(aes(x=abs, y=emmean), method="lm", se=FALSE)+
  scale_y_continuous(limits=c(-0.6,0.6))+
  scale_x_continuous(limits=c(-2.5,2.5))+
    ggtitle("Wellspring")
plot_W_fitness
```

![](README_files/figure-gfm/Fitness%20regression-2.png)<!-- -->

``` r
plot_2020_fitness <- plot_grid(plot_C_fitness, plot_W_fitness, labels="AUTO", align='h')
plot_2020_fitness
```

![](README_files/figure-gfm/Fitness%20regression-3.png)<!-- -->

``` r
save_plot("Figure_6.pdf", plot_2020_fitness, base_height=6, base_width=12)
```
