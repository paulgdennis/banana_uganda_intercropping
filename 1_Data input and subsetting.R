#############################################
# Georgina - Banana inter-cropping - Uganda #
#############################################

# Paul Dennis 2023

##############################
# Data input and preparation #
##############################

source('Functions.R')

### Load the OTU table with taxonomy into the environment - Change the file name and path to match yours
otu.tmp <- read.table("../Data/otu_with_tax_1750.csv",header=TRUE,sep=',',row.names = 1)

# Make Taxonomy table
taxonomy.all <- data.frame(OTU = row.names(otu.tmp),
                           Silva138_taxonomy = otu.tmp$taxonomy)

taxonomy.all <- data.frame(taxonomy.all,
                           separate(taxonomy.all, Silva138_taxonomy, sep = ";",
                                    c("Domain","Phylum", "Class",
                                      "Order", "Family", "Genus", "Species")))

taxonomy.all <- data.frame(taxonomy.all, 
                           d_p_           = paste(taxonomy.all$Domain,
                                                  taxonomy.all$Phylum, 
                                                  sep = ";"),
                           d_p_c_         = paste(taxonomy.all$Domain,
                                                  taxonomy.all$Phylum,
                                                  taxonomy.all$Class, 
                                                  sep = ";"),
                           d_p_c_o_       = paste(taxonomy.all$Domain,
                                                  taxonomy.all$Phylum,
                                                  taxonomy.all$Class,
                                                  taxonomy.all$Order,
                                                  sep = ";"),
                           d_p_c_o_f_     = paste(taxonomy.all$Domain,
                                                  taxonomy.all$Phylum,
                                                  taxonomy.all$Class,
                                                  taxonomy.all$Order,
                                                  taxonomy.all$Family,
                                                  sep = ";"),
                           d_p_c_o_f_g_   = paste(taxonomy.all$Domain,
                                                  taxonomy.all$Phylum,
                                                  taxonomy.all$Class,
                                                  taxonomy.all$Order,
                                                  taxonomy.all$Family,
                                                  taxonomy.all$Genus,
                                                  sep = ";"),
                           d_p_c_o_f_g_s_ = paste(taxonomy.all$Domain,
                                                  taxonomy.all$Phylum,
                                                  taxonomy.all$Class,
                                                  taxonomy.all$Order,
                                                  taxonomy.all$Family,
                                                  taxonomy.all$Genus,
                                                  taxonomy.all$Species,
                                                  sep = ";"))

row.names(taxonomy.all) <- taxonomy.all$OTU


### Create a sample (n) x OTU (p) table with relative abundances by transposing *otu.tmp without taxonomy column
otu <- as.data.frame(t(otu.tmp[,-98])/1750)

### Load environmental metadata into environment - Change the file name and path to match yours
env <- read.table("../Data/env.csv", header=TRUE, sep=',', row.names = 1)

### Convert pH.Group to a factor to explicitly define it as a categorical predictor variable. 
env$Treatment.full <- factor(env$Treatment.full,
                             levels = c("T0.Ndizi.Ctrl","T2.Ndiizi.rhizo","T5.Ndizi.rhizo",
                                        "T4.Bogoya.Ctrl","T1.Bogoya.rhizo","T6.Bogoya.rhizo","T3.Nakitembe.Ctrl"))
env <- env[order(env$Treatment.full),]
otu <- otu[row.names(env),]

### Check samples are in the same order
row.names(otu) == row.names(env)

# We'll make a color palette to represent the pH of each sample in plots
mypal.treatment <- colorRampPalette(c('#b10026','#005a32','#fc4e2a','#88419d','#41ab5d','#fed976','#d9f0a3'))

#################################
#################################
## Analyses of alpha diversity ##
#################################
#################################

#########################################################
# Univariate measures of central tenancy and dispersion #
#########################################################

# Get a summary for each level of Landuse
stby(env, env$Treatment.full, descr)

#############################################################################
# Univariate visualisations and Hypothesis tests for categorical predictors #
#############################################################################

# Sobs
### Determine whether alpha diversity differs between land uses using ANOVA and posthoc
mod1 <- aov(Sobs ~ Treatment.full, data = env)
summary(mod1)
cld(lsmeans(mod1, ~ Treatment.full), Letters = letters)

# Shannon
### Determine whether alpha diversity differs between land uses using ANOVA and posthoc
mod1 <- aov(Shan ~ Treatment.full, data = env)
summary(mod1)
cld(lsmeans(mod1, ~ Treatment.full), Letters = letters)

# Chao1
### Determine whether alpha diversity differs between land uses using ANOVA and posthoc
mod1 <- aov(Chao1 ~ Treatment.full, data = env)
summary(mod1)
cld(lsmeans(mod1, ~ Treatment.full), Letters = letters)

# Faith's PD
### Determine whether alpha diversity differs between land uses using ANOVA and posthoc
mod1 <- aov(PD ~ Treatment.full, data = env)
summary(mod1)
cld(lsmeans(mod1, ~ Treatment.full), Letters = letters)

### Visualise the data using a barplot with Standard Errors of the means as error bars   
bargraph.CI(Treatment.full, Sobs, data = env, col = mypal.treatment(7), main = "Sobs") 
bargraph.CI(Treatment.full, Shan, data = env, col = mypal.treatment(7), main = "Shannon") 
bargraph.CI(Treatment.full, Chao1, data = env, col = mypal.treatment(7), main = "Chao1") 
bargraph.CI(Treatment.full, PD, data = env, col = mypal.treatment(7), main = "PD") 

### Relationship with severity
anova(lm(Severity ~ Sobs, data = env))
anova(lm(Severity ~ Shan, data = env))
anova(lm(Severity ~ Chao1, data = env))
anova(lm(Severity ~ PD, data = env))

par(mfrow = c(2,2))
plot(Severity ~ Sobs, data = env)
plot(Severity ~ Shan, data = env)
plot(Severity ~ Chao1, data = env)
plot(Severity ~ PD, data = env)
abline(lm(Severity ~ Sobs, data = env))

### Visualise the data using a barplot with Standard Deviation as error bars   
#bargraph.CI(Treatment.full, Shan, data = env, col = mypal.treatment(7),
#            ci.fun= function(x) c(mean(x)-sd(x), mean(x) + sd(x))) 
 

##################################
### Analyses of Beta diversity ###
##################################

# Using square root transformed OTU relative abundances (aka Hellinger transformation)

################################################################################################################
# Determine whether the composition of bacterial communities is significantly influenced by pH using PERMANOVA #
################################################################################################################

# As a numeric predictor (pH)...
adonis2(sqrt(otu) ~ Treatment.full, data = env, method = 'euc')

adonis2(sqrt(otu) ~ Severity*Treatment.full, data = env, method = 'euc')


############################################################
# View community similar/dissimilarity in ordination space # 
############################################################

### First we'll make a distance-based PCA ordination object.
otu.pca <- rda(sqrt(otu) ) # PCA gives a euclidean projection of the input, which in dbPCA is transformed (e.g. Hellinger)
axis.percent(otu.pca) # Returns the proportion (%) of total inertia (compositional dissimilarity) explained on each axis 
otu.rda <- rda(sqrt(otu) ~ Severity, env )

# A very basic plot can be obtained as follows:
plot(otu.pca, type = 'n', scaling = 3)
points(otu.pca, dis='sites', pch = 21, bg = mypal.treatment(7)[env$Treatment.full], scaling = 3, cex = 2)
legend("bottomright", legend=unique(env$Treatment.full), pch=19,
       col=mypal.treatment(7)[unique((env$Treatment.full))])

### Here is a much more complex but informative plot: 

### Start 
# set the constant variables used below then plot the axes with %variation explained
ord = otu.pca
scaling.val = 2

plot(ord, 
     type='n', scaling=scaling.val, 
     xlab=paste("db-PC1 (",axis.percent(
       ord)[[1]],"%)",sep=""),
     ylab=paste("db-PC2 (",axis.percent(
       ord)[[2]],"%)",sep=""))

# Add ellipses representing the 95% confidence interval between treatments
ordiellipse(ord, env$Treatment.full, 
            col=mypal.treatment(7), 
            draw = 'polygon', alpha = 0.5, lty=0,  kind ='sd',scaling=scaling.val) # can also use kind = 'se' etc...

# Plot the OTUs in grey - we can't read most of them anyway (this can also be omitted)
points(ord, dis='sp', pch=4, col='grey', cex=0.6, scaling=scaling.val)

# Plot the Samples colored by pH.Group (for pH as a continuous variable we'd need 42 colors instead)
points(ord, dis='sites', pch=21, scaling = scaling.val,
       bg = addTrans(mypal.treatment(7)[factor(env$Treatment.full)], 180), cex = 2)

# The next block of code simply highlights which OTUs are most discriminating in red and adds their OTU_ID codes      
sd.val = 7 # This is the number of standard deviations from the mean position of OTUs on each axis
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),],
       pch=4,col='darkred',cex=0.6)	
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),],
       pch=4,col='darkred',cex=0.6) 
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),],
       pch=4,col='darkred',cex=0.6)
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),],
       pch=4,col='darkred',cex=0.6) 
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)

# This adds the legend
legend("bottomright",legend=unique(env$Treatment.full), pch=19,
       col=mypal.treatment(7)[unique(factor(env$Treatment.full))])

### END

plot(envfit(otu.pca, env$Severity))
#############################################################
# View the relative abundance of dominant OTUs in a heatmap #
#############################################################

# Select OTUs based on a relative abundance threshold 
hm.tmp <- otu[,which(apply(otu,2,max)>=0.01)]
dim(hm.tmp) # too many OTUs

# limit to OTUs that have a mean relative abundance of >=1% within treatment
otu.mean <- aggregate(otu,list(env$Treatment.full), mean)
row.names(otu.mean) <- otu.mean[,1]
otu.mean <- otu.mean[,-1]

hm.tmp <- otu.mean[,which(apply(otu.mean,2,max)>=0.01)]
dim(hm.tmp) # 28 OTUs 

# Get the taxonomy information for the OTUs in the heatmap
otus.for.heatmap.tmp <- taxonomy.all[colnames(hm.tmp),]
otus.for.heatmap <- otus.for.heatmap.tmp[order(otus.for.heatmap.tmp$Silva138_taxonomy),]

hm.mean <- hm.tmp[,row.names(otus.for.heatmap)]

hm.all <- otu[,row.names(otus.for.heatmap)]

# Make a full.id label by concatenated the OTU_ID (in square brackets) with the taxonomy
otus.for.heatmap$full.id <- paste("[", 
                                  otus.for.heatmap$OTU, 
                                  "] ", 
                                  otus.for.heatmap$Silva138_taxonomy, 
                                  sep='')

# Make a heatmap
mypal.bw <- colorRampPalette(c("White","Black"))
heatmap(t(sqrt(hm.all)), revC = TRUE, 
        Colv = NA, Rowv = NA, col=mypal.bw(180),
        labRow = otus.for.heatmap$full.id, labCol = env$Treatment.full, scale = 'none',margins = c(6,30))
        

for(i in otus.for.heatmap$OTU){
  print(i)
  print(anova(lm(otu[,i] ~ env$Treatment.full)))
  bargraph.CI(x.factor = env$Treatment.full, response = otu[,i], main = i)
}



