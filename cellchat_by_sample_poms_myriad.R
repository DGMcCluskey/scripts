library(Seurat)
library(zellkonverter)
library(CellChat)
library(future)
library(randomForest)
library(pROC)
library(caret)
library(pheatmap)

setwd("/myriadfs/home/rmhadgm/Scratch")

poms <- readH5AD("poms.annotated.dan.h5ad")

poms.seurat <- as.Seurat(poms)
rm(poms)


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

saveRDS(cellchat.list, file = "cellchat_list_by_sample.rds", compress = F)

#pull out interactions
interaction_data_list <- list()

# Loop over each CellChat object
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
            "Maternal", "Maternal", "Fetal", "Fetal")
  gestional_week = c("")
)



# Join metadata to feature matrix
rf_data <- left_join(metadata, feature_matrix, by = "sample")

# Remove sample column (keep group as label)
rf_data$sample <- NULL

set.seed(123)

# Random Forest classification
colnames(rf_data) <- make.names(colnames(rf_data))
rf_data$group <- as.factor(rf_data$group)
rf_model <- randomForest(group ~ ., data = rf_data, importance = TRUE, ntree = 1000)

# View important features (interactions)
importance_df <- as.data.frame(importance(rf_model))
importance_df$interaction <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]
head(importance_df, 20)

importance_noHLA <- importance_df[ !grepl("HLA", importance_df$interaction), ]

top_features <- head(importance_noHLA, 15)

ggplot(top_features, aes(x = reorder(interaction, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Interaction", y = "Importance (MeanDecreaseGini)", title = "Top Interactions by Random Forest")+
  theme_bw()+theme(axis.text = element_text(size = 14, colour = "black"))


interactions_plot <- pivot_longer(rf_data, names_to = "interaction", values_to = "prob", 2:6068)

interactions_plot <- filter(interactions_plot, interaction %in% top_features$interaction)

ggplot(interactions_plot, aes(x = group, y = prob, fill = group))+geom_boxplot()+theme_bw()+
  facet_wrap(~interaction, scales = "free")




aggregate(. ~ group, data = rf_data[, c("group", top_features)], FUN = mean)

set.seed(123)
# Set up cross-validation
cv_ctrl <- trainControl(method = "cv", number = 5)

# Train with CV
rf_cv_model <- train(
  group ~ ., data = rf_data,
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
conf_mat <- confusionMatrix(pred, rf_data$group)
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
train_idx <- createDataPartition(rf_data$group, p = 0.7, list = FALSE)

train_data <- rf_data[train_idx, ]
test_data  <- rf_data[-train_idx, ]

rf_model <- randomForest(group ~ ., data = train_data, importance = TRUE, ntree = 1000)

# Predict group for test samples
predictions <- predict(rf_model, newdata = test_data)

# Confusion matrix
conf_matrix <- confusionMatrix(predictions, test_data$group)
print(conf_matrix)

pheatmap(t(rf_data[, -1]), annotation_col = data.frame(Group = rf_data$group))
