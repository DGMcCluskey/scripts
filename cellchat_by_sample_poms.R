library(Seurat)
library(CellChat)
library(future)
library(dplyr)
library(tidyr)
library(randomForest)
library(caret)
library(zellkonverter)
library(ggplot2)
library(httpgd)
library(tidyr)
library(dplyr)


setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Tim_data")
getwd()
poms <- readH5AD("poms.annotated.dan.h5ad")

poms.seurat <- as.Seurat(poms)
rm(poms)



table(poms.seurat$patient_nice, poms.seurat$gestation_week)

# 1. Set identity to patient/sample
Idents(poms.seurat) <- "patient_nice"

# 2. (Optional) Remove any unwanted samples like "Fetal_9"
poms.seurat <- subset(poms.seurat, idents = "Fetal_9", invert = TRUE)

# 3. Split into list by sample
seurat.list <- SplitObject(poms.seurat, split.by = "patient_nice")
rm(poms.seurat)
# 4. Load and subset the CellChat DB
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)

# 5. Set up parallel processing
plan("multisession", workers = 4)
options(future.seed = TRUE)

# 6. Initialize output list
cellchat.list <- list()

# 7. Loop over each sample
for (sample_id in names(seurat.list)) {
  cat("Processing:", sample_id, "\n")
  seurat.obj <- seurat.list[[sample_id]]
  
  # Skip if sample is too small
  if (ncol(seurat.obj) < 100) {
    warning(paste("Skipping", sample_id, "— too few cells"))
    next
  }
  
  # Create CellChat object
  chat <- createCellChat(seurat.obj, group.by = "basic_labels", assay = "originalexp")
  chat@DB <- CellChatDB.use
  #drop unused levels
  chat@idents <- droplevels(chat@idents)
  # Core CellChat workflow
  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  chat <- computeCommunProb(chat, type = "triMean", population.size = TRUE)
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP")
  
  # Save intermediate CellChat object
  saveRDS(chat, file = paste0("cellchat_", sample_id, ".rds"))
  
  # Add to list
  cellchat.list[[sample_id]] <- chat
}

cellchat.list <- readRDS(file = "cellchat_list_by_sample.rds")

# Loop over each CellChat object
interaction_data_list <- list() 


for (sample_name in names(cellchat.list)) {
  chat_obj <- cellchat.list[[sample_name]]
  df <- subsetCommunication(chat_obj)
  
  if (!is.null(df)) {
    df$interaction_id <- paste(df$source, df$target, df$interaction_name, sep = "_")
    interaction_data_list[[sample_name]] <- df[, c("interaction_id", "prob")]
  }
}

# Combine all into a single wide matrix
feature_matrix <- bind_rows(interaction_data_list, .id = "sample") %>%
  pivot_wider(names_from = interaction_id, values_from = prob, values_fill = 0)

unique(feature_matrix[, 1])

#get metadata to add
metadata <- data.frame(
  sample = c("Fetal_6", "Maternal_10", "Maternal_2", "Fetal_12", "Fetal_2", "Fetal_1", "Maternal_6",
             "Maternal_4", "Fetal_13", "Fetal_8", "Fetal_10", "Maternal_7", "Maternal_12",
             "Maternal_1", "Maternal_8", "Fetal_7", "Fetal_4"),
  group = c("Fetal", "Maternal", "Maternal", "Fetal", "Fetal", "Fetal", "Maternal",
            "Maternal", "Fetal", "Fetal", "Fetal", "Maternal", "Maternal",
            "Maternal", "Maternal", "Fetal", "Fetal"),
  gestional_week = c(32.7,25.4,25.3,23.4,25.3,29.6,32.7,
                     35.9,25,34,25.4,24.1,23.4,
                     29.6,34,24.1,35.9),
  gestation_bin = c("moderate", "extreme", "extreme", "extreme", "extreme", "moderate", "moderate",
                    "moderate", "extreme", "moderate", "extreme", "extreme", "extreme",
                    "moderate", "moderate", "extreme", "moderate") 
  # using eliz's definitions for extremely, very and moderately preterm (moved the one "very" to moderate)
)

# Join metadata to feature matrix
rf_data <- left_join(metadata, feature_matrix, by = "sample")

# Remove sample column (keep group as label)
rf_data$sample <- NULL

set.seed(123)

#choose specific data - here only taking fetal samples to then look at time
rf_data <- filter(rf_data, group == "Fetal")
#remove unused meta columns, as we don't want those used as factors for the random forest
rf_data <- rf_data[, -c(1,2)]

# Random Forest classification
colnames(rf_data) <- make.names(colnames(rf_data))
#rf_data$group <- as.factor(rf_data$group)
rf_data$gestation_bin <- as.factor(rf_data$gestation_bin)
rf_model <- randomForest(gestation_bin ~ ., data = rf_data, importance = TRUE, ntree = 1000)

plot(rf_model$err.rate[,1], type = "l", ylab = "OOB Error", xlab = "Number of Trees")
plot(rf_model, main = "OOB Error Rate")
varImpPlot(rf_model, main = "Feature Importance")
#importance(rf_model)

# View important features (interactions)
importance_df <- as.data.frame(importance(rf_model))
importance_df$interaction <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$MeanDecreaseAccuracy), ]
head(importance_df, 20)

importance_noHLA <- importance_df[ !grepl("HLA", importance_df$interaction), ]

top_features <- head(importance_noHLA, 20)

p <- ggplot(top_features, aes(x = reorder(interaction, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy, fill = MeanDecreaseGini)) +
  geom_point(shape = 21, size = 6)+
  coord_flip() +
  labs(x = "Interaction", y = "Importance (MeanDecreaseAccuracy)", title = "Top interactions that distinguish\nmoderately and extremely preterm samples")+
  theme_bw()+theme(axis.text = element_text(size = 14, colour = "black"))+
  scale_fill_viridis_c()
p

#ggsave(plot = p, )

interactions_plot <- pivot_longer(rf_data, names_to = "interaction", values_to = "prob", 2:6068)

top_features <- head(importance_noHLA, 16)

interactions_plot <- filter(interactions_plot, interaction %in% top_features$interaction)

options(scipen = 999)

interactions_plot$gestation_bin <- factor(interactions_plot$gestation_bin, levels = c("extreme", "moderate"),
                                          labels = c("<28 gw", ">28 gw"))

ggplot(interactions_plot, aes(x = gestation_bin, y = prob, fill = gestation_bin))+
  geom_boxplot(outlier.shape = NA, alpha = 0.4, colour = "black")+
  geom_jitter(width = 0.1, shape = 21, colour = "black", size = 2)+
  theme_bw()+
  facet_wrap(~interaction, scales = "free_y", ncol = 4)+
  scale_fill_manual(values = c("darkorchid", "mediumseagreen"))+
  xlab(NULL)+
  theme(axis.text = element_text(size = 12, colour = "black"))

results <- interactions_plot %>%
  group_by(interaction) %>%
  summarise(
    p_value = wilcox.test(prob ~ gestation_bin)$p.value,
    .groups = "drop"
  ) %>%
  arrange(p_value)



set.seed(123)
# Set up cross-validation
cv_ctrl <- trainControl(method = "cv", number = 5)

# Train with CV
rf_cv_model <- train(
  gestation_bin ~ ., data = rf_data,
  method = "rf",
  trControl = cv_ctrl,
  importance = TRUE,
  ntree = 1000
)

# Accuracy, confusion matrix, etc.
print(rf_cv_model)

rf_cv_model$results       # Accuracy, Kappa
rf_cv_model$finalModel    # The random forest itself
varImp(rf_cv_model)   

# Predict on training data (CV is internal to caret)
pred <- predict(rf_cv_model, rf_data)

# Confusion matrix
conf_mat <- confusionMatrix(pred, rf_data$gestation_bin)
print(conf_mat)

# Create a confusion matrix table
cm_table <- as.data.frame(conf_mat$table)
colnames(cm_table) <- c("Reference", "Prediction", "Freq")

# Plot
ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix Heatmap")

#running model again but splitting training and test data
# Stratified sampling: maintains the group balance
train_idx <- createDataPartition(rf_data$gestation_bin, p = 0.7, list = FALSE)

train_data <- rf_data[train_idx, ]
test_data  <- rf_data[-train_idx, ]

rf_model <- randomForest(gestation_bin ~ ., data = train_data, importance = TRUE, ntree = 1000)

# Predict group for test samples
predictions <- predict(rf_model, newdata = test_data)

# Confusion matrix
conf_matrix <- confusionMatrix(predictions, test_data$gestation_bin)
print(conf_matrix)
library(pheatmap)
pheatmap(t(rf_data[, -1]), annotation_col = data.frame(Group = rf_data$gestation_bin))
