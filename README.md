Single versus mixed inoculations of bacteria on Lemna minor
================
Jason Laurich and Megan Frederickson
Sys.Date()

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
library(knitr) # render markdown
library(kableExtra) # exports tables to LaTeX format
```

### SECTION 2: Upload data

``` r
#####################Load 2020 data#######################
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

# Log abs count
df.CW$logabs<-log(df.CW$abs)

#Let's create separate dataframes for each population # Going to get the label 2020
df.C.2020<-subset(df.CW, df.CW$pop== "Churchill")
df.W.2020<-subset(df.CW, df.CW$pop== "Wellspring")

# Change factor order of bacteria to reflect their effects on duckweed growth.
#First, let's change the order of bacteria, ranked by effect. 
df.C.2020$bac<- factor(df.C.2020$bac, levels = c("Control","1","8","6","3","7","10","9","4","5","2","All10"))
# Assign names from sequencing
levels(df.C.2020$bac) <- c("Control", "Flavobacterium \nsuccinicans 1","Bosea \nmassiliensis","Aeromonas \nsalmonicida","Ohtaekwangia \nkoreensis","Flavobacterium \nsuccinicans 2","Falsiroseomonas \nstagni","Parasediminibacterium \npaludis","Arcicella sp.","Microbacterium \noxydans","Pseudomonas \nprotogens","All 10 bacteria")

df.W.2020$bac<- factor(df.W.2020$bac, levels = c("Control","2","8","1","9","6","3","5","4","7","10","All10"))
levels(df.W.2020$bac) <- c("Control","Sphingomonas \npituitosa 1","Flaviflagellibacter \ndeserti","Rhizobium \nrosettiformans","Rhizorhabdus \nwittichii 1","Rhizobium \ncapsici 1","Pseudomonas \nprotogens","Rhizorhabdus \nwittichii 2","Rhizobium \ncapsici 2","Sphingomonas \npituitosa 2","Fervidobacterium \nriparium","All 10 bacteria")

df2020 <- rbind(df.C.2020, df.W.2020)
df2020$pop <- trimws(df2020$pop)
```

### SECTION 3: Microbial growth models and figures

First, let’s fit a model combining all data (with bacterial strain
nested within population as a random effect) in which we test whether
the productivity of 10-strain communities differs from the productivity
of single strains. Again, this tests whether there is a
biodiversity-ecosystem function relationship, whereby the number of
strains in the microbiome (1 or 10) predicts microbial productivity.

``` r
#Group by whether 0, 1, or 10 strains were added
df2020$group <- as.factor(ifelse(df2020$bac == "All 10 bacteria", "10-strain community", ifelse(df2020$bac == "Control", "Control", "Single strain")))

#Make single strain the reference
df2020 = df2020 %>% mutate(group = relevel(group, "Single strain"))

#hist(log(df2020$abs))
df2020$concat <- paste0(df2020$pop, df2020$bac) #Concatenate pop and bac names, since there are different bacteria in each population
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
    ## group10-strain community        2.6116     0.6263  22.3734   4.170 0.000387 ***
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
df2020 = df2020 %>% mutate(bac = relevel(bac, "All 10 bacteria"))

#Churchill 2020
#Model 
lmer.abs.C<-lmer(abs~bac*plt + (1|edge) + (1|plate), data=subset(df2020, bac != "Control" & pop=="Churchill"))
summary(lmer.abs.C)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs ~ bac * plt + (1 | edge) + (1 | plate)
    ##    Data: subset(df2020, bac != "Control" & pop == "Churchill")
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
    ##                                         Estimate Std. Error       df t value
    ## (Intercept)                              3258.23     705.68    71.01   4.617
    ## bacFlavobacterium \nsuccinicans 1       -3156.33     956.05   201.78  -3.301
    ## bacBosea \nmassiliensis                 -3085.27    1035.73   201.64  -2.979
    ## bacAeromonas \nsalmonicida              -2570.84     935.65   200.92  -2.748
    ## bacOhtaekwangia \nkoreensis             -2787.08     976.94   201.32  -2.853
    ## bacFlavobacterium \nsuccinicans 2       -2713.76     934.61   201.05  -2.904
    ## bacFalsiroseomonas \nstagni             -3107.20     903.02   200.77  -3.441
    ## bacParasediminibacterium \npaludis       -683.37    1128.84   199.50  -0.605
    ## bacArcicella sp.                        -2407.99     934.97   201.82  -2.575
    ## bacMicrobacterium \noxydans             -3197.74     934.98   201.13  -3.420
    ## bacPseudomonas \nprotogens               4649.71    1000.20   199.20   4.649
    ## pltY                                     1828.70    1030.30   195.61   1.775
    ## bacFlavobacterium \nsuccinicans 1:pltY    908.95    1475.00   199.76   0.616
    ## bacBosea \nmassiliensis:pltY             -510.69    1435.68   200.69  -0.356
    ## bacAeromonas \nsalmonicida:pltY          -519.31    1331.92   197.86  -0.390
    ## bacOhtaekwangia \nkoreensis:pltY        -1466.68    1462.92   200.82  -1.003
    ## bacFlavobacterium \nsuccinicans 2:pltY    418.36    1407.16   196.32   0.297
    ## bacFalsiroseomonas \nstagni:pltY          399.04    1330.96   199.88   0.300
    ## bacParasediminibacterium \npaludis:pltY   -33.79    1562.85   197.54  -0.022
    ## bacArcicella sp.:pltY                      78.79    1374.16   197.73   0.057
    ## bacMicrobacterium \noxydans:pltY          322.86    1406.75   195.61   0.230
    ## bacPseudomonas \nprotogens:pltY          -411.63    1478.65   199.52  -0.278
    ##                                         Pr(>|t|)    
    ## (Intercept)                             1.69e-05 ***
    ## bacFlavobacterium \nsuccinicans 1       0.001138 ** 
    ## bacBosea \nmassiliensis                 0.003248 ** 
    ## bacAeromonas \nsalmonicida              0.006548 ** 
    ## bacOhtaekwangia \nkoreensis             0.004785 ** 
    ## bacFlavobacterium \nsuccinicans 2       0.004101 ** 
    ## bacFalsiroseomonas \nstagni             0.000705 ***
    ## bacParasediminibacterium \npaludis      0.545618    
    ## bacArcicella sp.                        0.010725 *  
    ## bacMicrobacterium \noxydans             0.000758 ***
    ## bacPseudomonas \nprotogens              6.07e-06 ***
    ## pltY                                    0.077466 .  
    ## bacFlavobacterium \nsuccinicans 1:pltY  0.538438    
    ## bacBosea \nmassiliensis:pltY            0.722429    
    ## bacAeromonas \nsalmonicida:pltY         0.697032    
    ## bacOhtaekwangia \nkoreensis:pltY        0.317273    
    ## bacFlavobacterium \nsuccinicans 2:pltY  0.766547    
    ## bacFalsiroseomonas \nstagni:pltY        0.764631    
    ## bacParasediminibacterium \npaludis:pltY 0.982774    
    ## bacArcicella sp.:pltY                   0.954337    
    ## bacMicrobacterium \noxydans:pltY        0.818715    
    ## bacPseudomonas \nprotogens:pltY         0.781010    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#plot(lmer.abs.C) # Not great, but easier/more interpretable to leave this on a raw scale for estimating additive expectation, rather than log-transforming response variable

#Additive and mean predictions
C_emmeans <- as.data.frame(emmeans(lmer.abs.C, specs = ~bac+plt))
C_micr_mean_noplt <- mean(C_emmeans[C_emmeans$plt == "N" & C_emmeans$bac != "All 10 bacteria", "emmean"])
C_micr_mean_plt <- mean(C_emmeans[C_emmeans$plt == "Y" & C_emmeans$bac != "All 10 bacteria", "emmean"])
C_micr_sum_noplt <- sum(C_emmeans[C_emmeans$plt == "N" & C_emmeans$bac != "All 10 bacteria", "emmean"])
C_micr_sum_plt <- sum(C_emmeans[C_emmeans$plt == "Y" & C_emmeans$bac != "All 10 bacteria", "emmean"])

#Wellspring 2020

#Model 
lmer.abs.W<-lmer(abs~bac*plt + (1|edge) + (1|plate), data=subset(df2020, bac != "Control" & pop=="Wellspring"))
summary(lmer.abs.W)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs ~ bac * plt + (1 | edge) + (1 | plate)
    ##    Data: subset(df2020, bac != "Control" & pop == "Wellspring")
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
    ##                                        Estimate Std. Error        df t value
    ## (Intercept)                            4290.429    642.371     7.826   6.679
    ## bacPseudomonas \nprotogens            -2760.081    846.801   183.242  -3.259
    ## bacSphingomonas \npituitosa 1         -3404.490    796.308   186.790  -4.275
    ## bacFlaviflagellibacter \ndeserti      -4127.509    793.843   186.568  -5.199
    ## bacRhizobium \nrosettiformans         -3307.629    761.587   179.633  -4.343
    ## bacRhizorhabdus \nwittichii 1         -4198.968    704.070   182.875  -5.964
    ## bacRhizobium \ncapsici 1              -4117.662    809.635   181.692  -5.086
    ## bacRhizorhabdus \nwittichii 2         -4098.765    820.393   186.515  -4.996
    ## bacRhizobium \ncapsici 2              -4150.789    788.669   185.392  -5.263
    ## bacSphingomonas \npituitosa 2         -3720.620    762.618   182.181  -4.879
    ## bacFervidobacterium \nriparium        -3339.743    765.632   181.802  -4.362
    ## pltY                                   2495.212    948.156   181.346   2.632
    ## bacPseudomonas \nprotogens:pltY        -238.306   1281.919   185.655  -0.186
    ## bacSphingomonas \npituitosa 1:pltY     -824.553   1253.536   186.876  -0.658
    ## bacFlaviflagellibacter \ndeserti:pltY -1627.915   1299.809   182.195  -1.252
    ## bacRhizobium \nrosettiformans:pltY      -23.689   1267.062   185.214  -0.019
    ## bacRhizorhabdus \nwittichii 1:pltY      799.189   1199.187   181.789   0.666
    ## bacRhizobium \ncapsici 1:pltY         -1514.928   1218.086   177.988  -1.244
    ## bacRhizorhabdus \nwittichii 2:pltY      176.455   1291.634   184.944   0.137
    ## bacRhizobium \ncapsici 2:pltY          2151.108   1217.382   180.957   1.767
    ## bacSphingomonas \npituitosa 2:pltY     3789.222   1221.788   182.488   3.101
    ## bacFervidobacterium \nriparium:pltY     974.289   1217.264   183.511   0.800
    ##                                       Pr(>|t|)    
    ## (Intercept)                           0.000172 ***
    ## bacPseudomonas \nprotogens            0.001331 ** 
    ## bacSphingomonas \npituitosa 1         3.04e-05 ***
    ## bacFlaviflagellibacter \ndeserti      5.22e-07 ***
    ## bacRhizobium \nrosettiformans         2.34e-05 ***
    ## bacRhizorhabdus \nwittichii 1         1.25e-08 ***
    ## bacRhizobium \ncapsici 1              9.07e-07 ***
    ## bacRhizorhabdus \nwittichii 2         1.34e-06 ***
    ## bacRhizobium \ncapsici 2              3.89e-07 ***
    ## bacSphingomonas \npituitosa 2         2.32e-06 ***
    ## bacFervidobacterium \nriparium        2.16e-05 ***
    ## pltY                                  0.009229 ** 
    ## bacPseudomonas \nprotogens:pltY       0.852728    
    ## bacSphingomonas \npituitosa 1:pltY    0.511488    
    ## bacFlaviflagellibacter \ndeserti:pltY 0.212020    
    ## bacRhizobium \nrosettiformans:pltY    0.985104    
    ## bacRhizorhabdus \nwittichii 1:pltY    0.505973    
    ## bacRhizobium \ncapsici 1:pltY         0.215247    
    ## bacRhizorhabdus \nwittichii 2:pltY    0.891485    
    ## bacRhizobium \ncapsici 2:pltY         0.078915 .  
    ## bacSphingomonas \npituitosa 2:pltY    0.002232 ** 
    ## bacFervidobacterium \nriparium:pltY   0.424519    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Additive and mean predictions
W_emmeans <- as.data.frame(emmeans(lmer.abs.W, specs = ~bac+plt))
W_micr_mean_noplt <- mean(W_emmeans[W_emmeans$plt == "N", "emmean"])
W_micr_mean_plt <- mean(W_emmeans[W_emmeans$plt == "Y", "emmean"])
W_micr_sum_noplt <- sum(W_emmeans[W_emmeans$plt == "N", "emmean"])
W_micr_sum_plt <- sum(W_emmeans[W_emmeans$plt == "Y", "emmean"])
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

# Generate Figure 2
Sum_2020_C_micr <- df.C.2020 %>% group_by(bac, plt) %>% summarize(n=n(), mean=mean(abs, na.rm=TRUE),
sd = sd(abs, na.rm=TRUE), sem = sd/sqrt(n))

#Add colors
Sum_2020_C_micr$color <- ifelse(Sum_2020_C_micr$bac == "All 10 bacteria", "10-strain community", ifelse(Sum_2020_C_micr$bac == "Control", "Control", "Single strain"))

#Churchill WITHOUT PLANTS
plot_C2020_micr_noplt <- ggplot(subset(Sum_2020_C_micr, bac != "Control" & plt == "N"), aes(y=mean,x=reorder(bac, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
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
   scale_x_discrete(labels=c("*Microbacterium oxydans*", "*Bosea massiliensis*",  "*Flavobacterium succinicans* 1", "*Falsiroseomonas stagni*", "*Ohtaekwangia koreensis*", "*Flavobacterium succinicans* 2", "*Aeromonas salmonicida*", "*Arcicella* sp.", "*Parasediminibacterium paludis*", "All 10 bacteria", "*Pseudomonas protogens*"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_C2020_micr_noplt_updated_v2 <-plot_C2020_micr_noplt_updated+facet_zoom(ylim = c(0, 10000))

#Churchill WITH PLANTS
plot_C2020_micr_plt <- ggplot(subset(Sum_2020_C_micr, bac != "Control" & plt == "Y"), aes(y=mean,x=reorder(bac, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
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
   scale_x_discrete(labels=c("*Ohtaekwangia koreensis*", "*Bosea massiliensis*", "*Aeromonas salmonicida*", "*Microbacterium oxydans*",  "*Arcicella* sp.", "*Falsiroseomonas stagni*", "*Flavobacterium succinicans* 2",  "*Flavobacterium succinicans* 1", "*Parasediminibacterium paludis*", "All 10 bacteria", "*Pseudomonas protogens*"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_C2020_micr_plt_updated_v2 <-plot_C2020_micr_plt_updated+facet_zoom(ylim = c(0, 10000))

# Wellspring
Sum_2020_W_micr <- df.W.2020 %>% group_by(bac, plt) %>% summarize(n=n(), mean=mean(abs, na.rm=TRUE),
sd = sd(abs, na.rm=TRUE), sem = sd/sqrt(n))

#Add colors
Sum_2020_W_micr$color <- ifelse(Sum_2020_W_micr$bac == "All 10 bacteria", "10-strain community", ifelse(Sum_2020_W_micr$bac == "Control", "Control", "Single strain"))

#Wellspring WITHOUT PLANTS
plot_W2020_micr_noplt <- ggplot(subset(Sum_2020_W_micr, bac != "Control" & plt == "N"), aes(y=mean,x=reorder(bac, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
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
   scale_x_discrete(labels=c("*Flaviflagellibacter deserti*", "*Rhizorhabdus wittichii* 2", "*Rhizorhabdus wittichii* 1", "*Rhizobium capsici* 2", "*Rhizobium capsici* 1",  "*Sphingomonas pituitosa* 2", "*Sphingomonas pituitosa* 1", "*Fervidobacterium riparium*", "*Rhizobium rosettiformans*", "*Pseudomonas protogens*", "All 10 bacteria"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_W2020_micr_noplt_updated_v2 <-plot_W2020_micr_noplt_updated+facet_zoom(ylim = c(0, 8000))

#Wellspring WITH PLANTS
plot_W2020_micr_plt <- ggplot(subset(Sum_2020_W_micr, bac != "Control" & plt == "Y"), aes(y=mean,x=reorder(bac, mean), color=color))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac, mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
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
   scale_x_discrete(labels=c("*Rhizobium capsici* 1", "*Flaviflagellibacter deserti*", "*Sphingomonas pituitosa* 1", "*Rhizorhabdus wittichii* 2", "*Rhizobium rosettiformans*", "*Rhizorhabdus wittichii* 1", "*Pseudomonas protogens*", "*Fervidobacterium riparium*", "*Rhizobium capsici* 2", "All 10 bacteria", "*Sphingomonas pituitosa* 2"))+
  theme(axis.text.x = ggtext::element_markdown())

plot_W2020_micr_plt_updated_v2 <-plot_W2020_micr_plt_updated+facet_zoom(ylim = c(0, 8000))

final_microbes <- plot_grid(plot_C2020_micr_noplt_updated_v2, plot_C2020_micr_plt_updated_v2,plot_W2020_micr_noplt_updated_v2, plot_W2020_micr_plt_updated_v2, nrow=2, labels="AUTO")
final_microbes
```

![](README_files/figure-gfm/Microbial%20growth%20figures-1.png)<!-- -->

``` r
save_plot("Figure_1.pdf", final_microbes, base_height=12, base_width=18)
```

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
C_alt <-ggplot(subset(Sum_2020_C_micr_wide, bac != "Control"))+geom_point(aes(x=mean_N, y=mean_Y, color=color))+geom_abline(aes(intercept=0, slope=1), linetype="dotted")+theme_cowplot()+xlab(expression("Microbial density without hosts (cells/µL)"))+ylab(expression("Microbial density with hosts (cells/µL)"))+
  geom_errorbar(aes(x=mean_N, ymax=mean_Y+sem_Y, ymin=mean_Y-sem_Y, color=color))+
    geom_errorbarh(aes(xmax=mean_N+sem_N, xmin=mean_N-sem_N, y=mean_Y, color=color))+
 theme(legend.title = element_blank())+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.55, 0.15))+
    ggtitle("Churchill")+
    geom_text_repel(data=subset(Sum_2020_C_micr_wide, bac != "Control" & bac!="All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac))),fontface = "italic", size=2.5)+
    geom_text_repel(data=subset(Sum_2020_C_micr_wide, bac == "All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac))), size=2.5)+
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
W_alt <-ggplot(subset(Sum_2020_W_micr_wide, bac != "Control"))+geom_point(aes(x=mean_N, y=mean_Y, color=color))+geom_abline(aes(intercept=0, slope=1), linetype="dotted")+theme_cowplot()+xlab(expression("Microbial density without hosts (cells/µL)"))+ylab(expression("Microbial density with hosts (cells/µL)"))+
 theme(legend.title = element_blank())+
    geom_errorbar(aes(x=mean_N, ymax=mean_Y+sem_Y, ymin=mean_Y-sem_Y, color=color))+
    geom_errorbarh(aes(xmax=mean_N+sem_N, xmin=mean_N-sem_N, y=mean_Y, color=color))+
  scale_color_manual(values=c(dot_colors[1], dot_colors[3]))+
  theme(legend.position = c(0.55, 0.15))+
    ggtitle("Wellspring")+
  geom_text_repel(data=subset(Sum_2020_W_micr_wide, bac != "Control" & bac!="All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac))),fontface = "italic", size=2.5)+
    geom_text_repel(data=subset(Sum_2020_W_micr_wide, bac == "All 10 bacteria"), aes(x=mean_N, y=mean_Y, label=trimws(gsub("[\r\n]", "", bac))), size=2.5)+
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
save_plot("Figure_2.pdf", alt, base_height=6, base_width=12)
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
lm_both_frd <- lmer(frd1~frd0+group+(1|edge)+(1|plate)+(1|pop/bac), data=df2020)
summary(lm_both_frd)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: frd1 ~ frd0 + group + (1 | edge) + (1 | plate) + (1 | pop/bac)
    ##    Data: df2020
    ## 
    ## REML criterion at convergence: 2377.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0495 -0.6193  0.0055  0.6886  3.4129 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  plate    (Intercept) 0.44649  0.6682  
    ##  bac:pop  (Intercept) 0.03772  0.1942  
    ##  pop      (Intercept) 0.48787  0.6985  
    ##  edge     (Intercept) 0.26278  0.5126  
    ##  Residual             8.31245  2.8831  
    ## Number of obs: 475, groups:  plate, 39; bac:pop, 24; pop, 2; edge, 2
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)                4.8002     0.9011   6.9396   5.327  0.00112 ** 
    ## frd0                       2.0736     0.2395 464.1434   8.657  < 2e-16 ***
    ## groupSingle strain         0.4275     0.5077  19.6950   0.842  0.40983    
    ## group10-strain community   1.7780     0.6920  20.5294   2.570  0.01806 *  
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
    ## (Intercept) 28.3781  1  9.979e-08 ***
    ## frd0        74.9396  1  < 2.2e-16 ***
    ## group        8.0413  2    0.01794 *  
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
lm_both_pix <- lmer(pix1~scale(pix0)+group+(1|edge)+(1|plate)+(1|pop/bac), data=df2020)
summary(lm_both_pix)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pix1 ~ scale(pix0) + group + (1 | edge) + (1 | plate) + (1 |  
    ##     pop/bac)
    ##    Data: df2020
    ## 
    ## REML criterion at convergence: 9985.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.76518 -0.54633  0.05158  0.62057  2.87974 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  plate    (Intercept) 10596583 3255    
    ##  bac:pop  (Intercept)   379421  616    
    ##  pop      (Intercept) 31211567 5587    
    ##  edge     (Intercept)  8116286 2849    
    ##  Residual             82393701 9077    
    ## Number of obs: 475, groups:  plate, 39; bac:pop, 24; pop, 2; edge, 2
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
    ## pix1 ~ scale(pix0) + group + (1 | edge) + (1 | plate) + (1 | bac:pop) + (1 | pop)
    ##               npar  logLik   AIC     LRT Df Pr(>Chisq)    
    ## <none>           9 -4992.8 10004                          
    ## (1 | edge)       8 -5000.8 10018 16.0662  1  6.117e-05 ***
    ## (1 | plate)      8 -5002.3 10021 19.0147  1  1.297e-05 ***
    ## (1 | bac:pop)    8 -4992.8 10002  0.0683  1     0.7938    
    ## (1 | pop)        8 -5006.9 10030 28.3131  1  1.032e-07 ***
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
df2020 = df2020 %>% mutate(bac = relevel(bac, "Control"))

#Churchill
C_pixels <- lmer(pix1~pix0+bac+(1|edge)+(1|plate), data=subset(df2020, pop == "Churchill"))
summary(C_pixels)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pix1 ~ pix0 + bac + (1 | edge) + (1 | plate)
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
    ##                                     Estimate Std. Error        df t value
    ## (Intercept)                        3878.8193  3030.8194    4.8295   1.280
    ## pix0                                  4.9076     0.5278  220.0131   9.298
    ## bacAll 10 bacteria                 7874.6934  2406.2248  214.3893   3.273
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
    ##                                    Pr(>|t|)    
    ## (Intercept)                         0.25865    
    ## pix0                                < 2e-16 ***
    ## bacAll 10 bacteria                  0.00124 ** 
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
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#plot(C_pixels) #Looks ok

#Calculate additive and mean predictions 
C_pix_emmeans <- as.data.frame(emmeans(C_pixels, specs=~bac))
pred.C.2020.avg.emmeans <- mean(as.numeric(C_pix_emmeans[3:12, "emmean"])) 
pred.C.2020.sum.emmeans <- sum(as.numeric(C_pix_emmeans[3:12, "emmean"]-as.numeric(C_pix_emmeans[1, "emmean"])))+as.numeric(C_pix_emmeans[1, "emmean"]) #Add the difference between the mean for each single strain and the control, and then add the control

#Wellspring
W_pixels <- lmer(pix1~pix0+bac+(1|edge)+(1|plate), data=subset(df2020, pop == "Wellspring"))
summary(W_pixels)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pix1 ~ pix0 + bac + (1 | edge) + (1 | plate)
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
    ##                                    Estimate Std. Error         df t value
    ## (Intercept)                       5044.0552  3623.6227    10.5412   1.392
    ## pix0                                 7.2913     0.6512   223.2800  11.196
    ## bacAll 10 bacteria                4497.5803  3521.5576   222.8469   1.277
    ## bacPseudomonas \nprotogens        4856.7650  3527.2609   223.3201   1.377
    ## bacSphingomonas \npituitosa 1     1484.5091  3560.5816   220.9672   0.417
    ## bacFlaviflagellibacter \ndeserti -1077.4527  3553.0940   221.5627  -0.303
    ## bacRhizobium \nrosettiformans      659.1365  3497.7890   219.4910   0.188
    ## bacRhizorhabdus \nwittichii 1     -581.6311  3490.8712   217.3165  -0.167
    ## bacRhizobium \ncapsici 1          1307.6764  3518.4586   223.0462   0.372
    ## bacRhizorhabdus \nwittichii 2     3747.2298  3502.0644   220.2646   1.070
    ## bacRhizobium \ncapsici 2          4053.7219  3500.9568   219.8402   1.158
    ## bacSphingomonas \npituitosa 2     5636.4924  3502.5996   220.2029   1.609
    ## bacFervidobacterium \nriparium    1449.3299  3580.4124   224.0491   0.405
    ##                                  Pr(>|t|)    
    ## (Intercept)                         0.193    
    ## pix0                               <2e-16 ***
    ## bacAll 10 bacteria                  0.203    
    ## bacPseudomonas \nprotogens          0.170    
    ## bacSphingomonas \npituitosa 1       0.677    
    ## bacFlaviflagellibacter \ndeserti    0.762    
    ## bacRhizobium \nrosettiformans       0.851    
    ## bacRhizorhabdus \nwittichii 1       0.868    
    ## bacRhizobium \ncapsici 1            0.710    
    ## bacRhizorhabdus \nwittichii 2       0.286    
    ## bacRhizobium \ncapsici 2            0.248    
    ## bacSphingomonas \npituitosa 2       0.109    
    ## bacFervidobacterium \nriparium      0.686    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling

``` r
#plot(W_pixels) #Looks ok

#Calculate additive and mean predictions 
#emmeans
W_pix_emmeans <- as.data.frame(emmeans(W_pixels, specs=~bac))
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

Sum_2020 <- df2020 %>% group_by(pop, bac, group) %>% summarize(n=n(), mean=mean(pix1, na.rm=TRUE), sd=sd(pix1, na.rm=TRUE), sem=sd/sqrt(n))

plot_C2020 <- ggplot(subset(Sum_2020, pop == "Churchill"), aes(y=mean,x=reorder(bac, mean), color=group))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac,mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
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

plot_C2020_updated <- plot_C2020+scale_x_discrete(labels=c("*Flavobacterium succinicans* 1", "*Bosea massiliensis*", "Control", "*Flavobacterium succinicans* 2", "*Aeromonas salmonicida*", "*Ohtaekwangia koreensis*", "*Falsiroseomonas stagni*", "*Parasediminibacterium paludis*", "*Arcicella* sp.", "*Microbacterium oxydans*", "All 10 bacteria",  "*Pseudomonas protogens*"))+
  theme(axis.text.x = ggtext::element_markdown())
plot_C2020_updated
```

![](README_files/figure-gfm/Duckweed%20growth%20figures-1.png)<!-- -->

``` r
#Wellspring 2020
plot_W2020 <- ggplot(subset(Sum_2020, pop == "Wellspring"), aes(y=mean,x=reorder(bac, mean), color=group))+
  geom_point(size=pt_size)+
  geom_errorbar(aes(x=reorder(bac,mean), ymin=mean-sem,ymax=mean+sem),width=0.5)+
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

plot_W2020_updated <- plot_W2020+scale_x_discrete(labels=c("*Sphingomonas pituitosa* 1", "Control",  "*Flaviflagellibacter deserti*", "*Rhizobium rosettiformans*", "*Rhizobium capsici* 1", "*Rhizorhabdus wittichii* 1", "*Pseudomonas protogens*", "*Rhizorhabdus wittichii* 2", "*Rhizobium capsici* 2", "*Sphingomonas pituitosa* 2", "All 10 bacteria", "*Fervidobacterium riparium*"))+
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
save_plot("Figure_3.pdf", plot_2020, base_height=8, base_width=12)
```

### SECTION 5: Fitness regression

Fit models (with full random effect structure) from which we can extract
‘genotype’ (i.e. bacterial inocula) means for their (1) effects on
duckweed growth, and (2) microbial density when grown with plants. Then,
we will use scaled data so that we can center fitness proxy data around
0.

``` r
#Churchill
df.C.2020 <- subset(df2020, bac != "Control" & plt == "Y" & pop == "Churchill")
df.C.2020$pix1_scaled <- scale(df.C.2020$pix1, center = TRUE, scale = TRUE)
df.C.2020$pix0_scaled <- scale(df.C.2020$pix0, center = TRUE, scale = TRUE)
df.C.2020$abs_scaled <- scale(df.C.2020$abs, center=TRUE, scale=TRUE)

# Create models to extract strain/genotype means
modC <-lmer(pix1_scaled ~ bac + pix0_scaled + (1|edge) +(1|plate), data=df.C.2020)
emC <- as.data.frame(emmeans(modC, "bac", var="pix1_scaled"))
emC$CI<-emC$SE*qnorm(0.975)

modC.abs<-lmer(abs_scaled ~ bac + (1|edge) +(1|plate), data=df.C.2020)
em.C.abs<-as.data.frame(emmeans(modC.abs, "bac", var="abs"))

emC$color <- ifelse(emC$bac == "All 10 bacteria", "10-strain community", "Single strain")

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
  geom_text_repel(data=subset(emC, bac != "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac))),fontface = "italic", size=2.5)+
      geom_text_repel(data=subset(emC, bac == "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac))), size=2.5)+
  geom_smooth(aes(x=abs, y=emmean), method="lm", se=FALSE)+
    scale_y_continuous(limits=c(-0.6,0.6))+
  scale_x_continuous(limits=c(-2,2.5))+
      ggtitle("Churchill")
plot_C_fitness
```

![](README_files/figure-gfm/Fitness%20regression-1.png)<!-- -->

``` r
#Wellspring
df.W.2020 <- subset(df2020, bac != "Control" & plt == "Y" & pop == "Wellspring")
df.W.2020$pix1_scaled <- scale(df.W.2020$pix1, center = TRUE, scale = TRUE)
df.W.2020$pix0_scaled <- scale(df.W.2020$pix0, center = TRUE, scale = TRUE)
df.W.2020$abs_scaled <- scale(df.W.2020$abs, center=TRUE, scale=TRUE)

# Create models to extract strain/genotype means
modW <-lmer(pix1_scaled ~ bac + pix0_scaled + (1|edge) +(1|plate), data=df.W.2020)
emW <- as.data.frame(emmeans(modW, "bac", var="pix1_scaled"))
emW$CI<-emW$SE*qnorm(0.975)

modW.abs<-lmer(abs_scaled ~ bac + (1|edge) +(1|plate), data=df.W.2020)
em.W.abs<-as.data.frame(emmeans(modW.abs, "bac", var="abs"))

emW$color <- ifelse(emW$bac == "All 10 bacteria", "10-strain community", "Single strain")

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
  geom_text_repel(data=subset(emW, bac != "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac))),fontface = "italic", size=2.5)+
      geom_text_repel(data=subset(emW, bac == "All 10 bacteria"), aes(x=abs, y=emmean, label=trimws(gsub("[\r\n]", "", bac))), size=2.5)+
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
save_plot("Figure_4.pdf", plot_2020_fitness, base_height=6, base_width=12)
```
