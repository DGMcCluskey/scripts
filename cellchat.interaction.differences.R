library(CellChat)
library(circlize)
library(tidyr)
library(dplyr)
library(pheatmap)
#code for Claire Smith (through Marko)

#pull out the underlying interactions for each CellChat object
your.interactions.1 <- subsetCommunication(your.cellchat.1, slot.name = "net")
your.interactions.2 <- subsetCommunication(your.cellchat.2, slot.name = "net")

#getting number of sent and received interactions for each pathway, per celltype
group1 <- your.interactions.1 %>% group_by(pathway_name, source, target) %>% summarise(number = n())
group2 <- your.interactions.2 %>% group_by(pathway_name, source, target) %>% summarise(number = n())

#renaming columns ready for merge
names(group1)[4] <- "number1"
names(group2)[4] <- "number2"

#merging datagrame into one
combined <- merge(group1, group2, all = T)
#chaning NA values to 0 (NAs are here if interaction wasn't present in one dataset)
combined[is.na(combined)] <- 0
#difference in interactions between groups
combined$difference <- combined$number1-combined$number2

#pull out pathway you want to inspect
your.pathway <- filter(combined, pathway_name %in% "COLLAGEN")
#keep only sender/target and difference columns
your.pathway <- your.pathway[, c(2,3,6)]

#making colour scale that will reflect change in interaction number
col_fun = colorRamp2(range(your.pathway$difference), c("royalblue", "red"), transparency = 0.5)
#plotting chord plot
chordDiagram(your.pathway, grid.col = "grey50", col = col_fun, 
             link.zindex = rank(selplg[[3]]))

#alternative visualisation methods
#slot.name determines whether you will pull out pathway level or interaction level interactions
#netP is pathway, net gets individual interactions
your.interactions.1 <- subsetCommunication(your.cellchat.1, slot.name = "net")
#counting the number of source and target interactions for each pathway and celltype
your.interactions.source.1 <- your.interactions.1 %>% group_by(pathway_name, source) %>% summarise(number = n())
your.interactions.target.1 <- your.interactions.1 %>% group_by(pathway_name, target) %>% summarise(number = n())

#do the same for the other group/condition
your.interactions.2 <- subsetCommunication(Your.cellchat.2, slot.name = "net")
your.interactions.source.2 <- your.interactions.2 %>% group_by(pathway_name, source) %>% summarise(number = n())
your.interactions.target.2 <- your.interactions.2 %>% group_by(pathway_name, target) %>% summarise(number = n())

#renaming columns
names(your.interactions.source.1)[3] <- "group1.source.number"
names(your.interactions.source.2)[3] <- "group2.source.number"
names(your.interactions.target.1)[3] <- "group1.target.number"
names(your.interactions.target.2)[3] <- "group2.target.number"


#merging both source dataframes into one
source.diff <- merge(your.interactions.source.1, your.interactions.source.2, all = T)
#chaning NA values to 0 (NAs are here if interaction wasn't present in one dataset)
source.diff[is.na(source.diff)] <- 0
#difference in interactions between groups
source.diff$source.difference <- source.diff$group1.source.number - source.diff$group2.source.number

#do the above for target interactions
target.diff <- merge(your.interactions.target.1, your.interactions.target.2, all = T)
target.diff[is.na(target.diff)] <- 0
target.diff$target.difference <- target.diff$group1.target.number - target.diff$group2.target.number

#making common celltype column name
names(source.diff)[2] <- "celltype"
names(target.diff)[2] <- "celltype"

#combining both dataframes into one for plots
interactions.differences <- merge(source.diff, target.diff, all = T)
#chanign NAs to 0
interactions.differences[is.na(interactions.differences)] <- 0

#at this point, you have the differences in interactions of each celltype as EITHER senders or targets
#lets make a column of the total differences in interactions regardless of being a sender or target of the interaction
#bare in mind this might make differences smaller as they might be an increased sender but a even more of a reduced target
interactions.differences$overall.difference <- interactions.differences$source.difference+interactions.differences$target.difference


#if you only want to pull out a specific pathway of interest for plotting
pathways.of.interest <- filter(interactions.differences, pathway_name %in% c("SELPLG", "IL6", "IL4","MIF"))

#or plot using all data
#plotting barplot
#x variable can either be source, target or overall difference
ggplot(pathways.of.interest, aes(x = overall.difference, y = pathway_name, fill = celltype))+
  geom_col(position = position_dodge(0.5), width = 0.5, colour = "black")+
  theme_bw()+geom_vline(xintercept = 0, linetype = "dashed", colour = "black")+
  xlab("Difference in number of interaction number group1 vs group2")+ylab("Pathway")

#plotting heatmap
heatmap.input <- pathways.of.interest[, c(1,2,9)]
heatmap.input <- pivot_wider(heatmap.input, names_from = "celltype", values_from = "overall.difference")
heatmap.input[is.na(heatmap.input)] <- 0
heatmap.input <- as.matrix(heatmap.input)
rownames(heatmap.input) <- heatmap.input[, 1]
heatmap.input <- heatmap.input[, -1]
class(heatmap.input) <- "numeric"

heatmap.input
pheatmap::pheatmap(heatmap.input, scale = "column", color = colorRampPalette(brewer.pal(10, "YlGnBu"))(50), 
                   border_color = "black",
                   cluster_cols = F, cellheight = 20, cellwidth = 20)
