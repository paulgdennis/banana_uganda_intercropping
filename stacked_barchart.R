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



corrplot(cor(otu_p_c[which(apply(otu_p_c, 2, max)>=0.01)]), diag = F, type = 'upper')

# Heatmap
mypal <- colorRampPalette(c("White", "Black"))

heatmap(as.matrix(t(sqrt(otu_p_c.mean))), 
        Rowv = NA, Colv = NA, scale = "none", revC = TRUE, col=mypal(256))

# Pie chart
pie(sort(colSums(otu_p_c),decreasing = TRUE), 
    clockwise = TRUE, main='Relative frequency - Phylum-Class')





