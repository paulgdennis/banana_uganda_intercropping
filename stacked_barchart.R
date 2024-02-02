### Stacked bargraph
library(corrplot)

## Relative frequencies of Phylum-Class 
sqrt(otu)
# Summarise frequencies 
otu_p_c <- aggregate(t(otu),
                     list(taxonomy.all[colnames(otu),]$d_p_c_), sum)
row.names(otu_p_c) <- otu_p_c[,1] 
otu_p_c <- as.data.frame(t(otu_p_c[,-1]))
otu_p_c
dim(otu_p_c)

otu_p_c.mean <- aggregate((otu_p_c), list(env$Treatment.full), mean)
row.names(otu_p_c.mean) <- otu_p_c.mean[,1] 
otu_p_c.mean <- as.data.frame((otu_p_c.mean[,-1]))
dim(otu_p_c.mean)

dim(otu_p_c.mean[,which(apply(otu_p_c.mean,2,max)>=0.005)]*100)
# Stacked bargraph

T0.Ndizi.Ctrl T2.Ndiizi.rhizo T5.Ndizi.rhizo T4.Bogoya.Ctrl T1.Bogoya.rhizo T6.Bogoya.rhizo T3.Nakitembe.Ctrl
data.frame(t(otu_p_c.mean[,which(apply(otu_p_c.mean,2,max)>=0.005)]*100)[,1]
           
rbind(t(otu_p_c.mean[,which(apply(otu_p_c.mean,2,max)>=0.005)]*100)[,1],
      t(otu_p_c.mean[,which(apply(otu_p_c.mean,2,max)>=0.005)]*100)[,2])

x <- data.frame(Group = names(otu_p_c.mean[,which(apply(otu_p_c.mean,2,max)>=0.005)]),
                Relative_abundance = 100*colSums(otu_p_c.mean))

df <- data.frame(x, separate(x, Group, sep = ";", 
                             c("Domain","Phylum","Class")))
df
ggnested(df, aes(x = Phylum, y = Relative_abundance, main_group = Phylum, sub_group = Class)) +
  geom_bar(stat = "identity", aes(x = Class))



#corrplot(cor(otu_p_c[which(apply(otu_p_c, 2, max)>=0.01)]), diag = F, type = 'upper')

# Heatmap
mypal <- colorRampPalette(c("White", "Black"))

heatmap(as.matrix(t(sqrt(otu_p_c))), 
        Rowv = NA, Colv = NA, scale = "none", revC = TRUE, col=mypal(256))

# Pie chart
pie(sort(colSums(otu_p_c),decreasing = TRUE), 
    clockwise = TRUE, main='Relative frequency - Phylum-Class')


## Alpha diversity of Supergroup;Divisions 

# Calculate Shannon
shan.all.7500 <- diversity(otu.all.7500,"shannon")

# Calculate Sobs
sobs.all.7500 <- rowSums(decostand(otu.all.7500, method='pa'))

sg_d_sd_.sobs.all.7500 <- aggregate(t(decostand(otu.all.7500, method='pa')),
                                    list(taxonomy.all[colnames(otu.all.7500),]$sg_d_sd_), sum)
row.names(sg_d_sd_.sobs.all.7500) <- sg_d_sd_.sobs.all.7500[,1] 
sg_d_sd_.sobs.all.7500 <- as.data.frame(t(sg_d_sd_.sobs.all.7500[,-1]))
sg_d_sd_.sobs.all.7500

# Stacked bargraph of Sobs
df <- data.frame(df, Richness = colSums(sg_d_sd_.sobs.all.7500))

ggnested(df, aes(x = Supergroup, y = Richness, main_group = Supergroup, sub_group = Subdivision)) +
  geom_bar(stat = "identity", aes(x = Supergroup))

# Pie chart
pie(sort(colSums(sg_d_sd_.sobs.all.7500),decreasing = TRUE), 
    clockwise = TRUE, main='Sobs - Supergroups')

## Gradient analyses

plot(sobs.all.7500 ~ env.all$Lat)
plot(shan.all.7500 ~ env.all$Lat)

plot(sobs.all.7500 ~ env.all$Temp)
plot(shan.all.7500 ~ env.all$Temp)








# Lollypop plots


# Libraries
library(tidyverse)
library(hrbrthemes)
library(kableExtra)
options(knitr.table.format = "html")
library(patchwork)

df %>%
  filter(!is.na(Richness)) %>%
  arrange(Richness) %>%
  tail(10) %>%
  mutate(Group=factor(Group, Group)) %>%
  ggplot( aes(x=Group, y=Richness) ) +
  geom_segment( aes(x=Group ,xend=Group, y=0, yend=Richness), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") 

df %>%
  filter(!is.na(Relative_abundance)) %>%
  arrange(Relative_abundance) %>%
  tail(10) %>%
  mutate(Group=factor(Group, Group)) %>%
  ggplot( aes(x=Group, y=Relative_abundance) ) +
  geom_segment( aes(x=Group ,xend=Group, y=0, yend=Relative_abundance), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") 
