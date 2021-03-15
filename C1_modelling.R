# Get data splits (train and test set) for the location grid cells.
# Preprocess the data (removing NAs etc.).
# Train the model, choosing from logistic regression, random forest and support vector machine.
# Use the final model to create an incident likelihood map.
# Measure model performance on the test set.
# Explore the relationship of the features with the dependent variable (incident yes/no).


# GIS libraries
library(raster)
# modelling libraries
library(tidyverse)
library(tidymodels)
library(caret)
# visualization libraries
library(corrplot)
library(ggplot2)
library(plotROC)


#------------------------------------------------------------
# setup
#------------------------------------------------------------

# choose a scenario name
sc_name <- 'SC_1month' # other scenarios: 'SC_6month', 'SC_3month', 'SC_alldata'

model_type <- 'rf' # other model types: 'log', 'svm'

inputs <- read.csv(paste0('output/location_features_', sc_name, '.csv'))

out <- stack(paste0('output/rasters_', sc_name, '.tif'))


#------------------------------------------------------------
# look into input locations and their corresponding features
#------------------------------------------------------------

summary(inputs)
str(inputs)
table(inputs$train)

# correlation check for predictors
corrplot(cor(inputs %>% 
               dplyr::select(- c(train, test, cell_id, train_label, test_label)) %>% 
               filter(! is.na(elevation_m) & ! is.na(slope_deg) & ! is.na(ndvi_mean)))) 


#------------------------------------------------------------
# set divisions
#------------------------------------------------------------

# the training and test set contain label 0/1 for absence/presence of an indicent
# and label 2 represents the remaining locations in the area of interest


train <- inputs %>% 
  filter(train != 2) %>% 
  mutate(train = factor(train, levels = c(0,1), labels = c("Unharmed", "Incident")))
dim(train)
summary(train)
table(train$train)


test <- inputs %>% 
  filter(test != 2) %>% 
  mutate(test = factor(test, levels = c(0,1), labels = c("Unharmed", "Incident")))
dim(test)
summary(test)
table(test$test)


#------------------------------------------------------------
# preprocessing
#------------------------------------------------------------

recipe <-
  train %>%
  recipe(train ~ .) %>%
  step_rm(test, train_label, test_label) %>% 
  step_naomit(all_predictors()) %>% 
  step_corr(all_numeric(), -all_outcomes(), threshold = 0.9) %>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  prep()
# step_normalize(all_numeric(), -all_outcomes()) %>% # grootte coefficenten vergelijken

baked_train <- bake(recipe, train)
dim(baked_train)
summary(baked_train)
table(baked_train$train)


#------------------------------------------------------------
# training the logistic regression model
#------------------------------------------------------------

if (model_type == 'log'){
  mod_basic <- logistic_reg(mode = "classification") %>%
    set_engine("glm") %>%
    fit(train ~ . -cell_id , data = baked_train) # -dtown -dwater -droad
  tidy(mod_basic) 
  
  
  # step() only accepts glm() model, it gives the same values as the tidymodels model
  mod_basic2 <- glm(train ~ . -cell_id, data = baked_train, family = 'binomial')
  summary(mod_basic2)
  
  # stepwise feature selection to determine best predictors
  stats::step(mod_basic2) # backwards feature selection
  
  
  # mod_final <- glm(formula = train ~ elevation_m + ndvi_mean + ndvi_max, 
  #                  family = "binomial", data = baked_train)
  # summary(mod_final)
  mod_final <- logistic_reg(mode = "classification") %>%
    set_engine("glm") %>%
    fit(train ~ dwater + droad, data = baked_train)
  tidy(mod_final)
  mod_final
}

if (model_type == 'rf'){
  mod_basic <- rand_forest(mode = "classification") %>%
    set_engine("randomForest") %>%
    fit(train ~ . -cell_id, data = baked_train)
  mod_basic
  
  # step() only accepts glm() model
  
  mod_final <- rand_forest(mode = "classification", 
                           mtry=5,
                           trees=100,
                           min_n=5) %>%
    set_engine("randomForest") %>%
    fit(train ~ . -cell_id, data = baked_train)
  mod_final
}

if (model_type == 'svm'){
  mod_basic <- svm_rbf(mode = "classification") %>%
    set_engine("kernlab") %>%
    fit(train ~ . -cell_id, data = baked_train)
  mod_basic
  
  # step() only accepts glm() model
  
  mod_final <- svm_rbf(mode = "classification", 
                       cost = 10, 
                       rbf_sigma = 0.035) %>%
    set_engine("kernlab") %>%
    fit(train ~ . -cell_id, data = baked_train)
  mod_final
}


# # cross-validation
# # ----------
# set.seed(555)
# l_cv <- vfold_cv(baked_train, v = 10, strata = "train")
# 
# mod_glm <-map2_df(.x = l_cv$splits,
#                   .y = l_cv$id,
#                   function (split = .x, fold = .y) 
#                   {
#                     # Split the data into analysis and assessment tables
#                     df_analysis <- analysis(split)
#                     df_assessment <- assessment(split)
#                       
#                     # Build the model
#                     mod <-
#                       logistic_reg(mode = "classification") %>%
#                       set_engine("glm") %>%
#                       fit(train ~ ., data = df_analysis)
#                       
#                     # Summarise Predictions
#                     table <- 
#                       tibble(fold = fold,
#                               truth = df_assessment$train,
#                               .pred_Incident = 
#                                 predict(mod, 
#                                         new_data = df_assessment, 
#                                         type = "prob")[[".pred_Incident"]],
#                               .pred_Unharmed = 
#                                 predict(mod, 
#                                         new_data = df_assessment, 
#                                         type = "prob")[[".pred_Unharmed"]],
#                               .pred_Class = 
#                                 predict(mod, 
#                                         new_data = df_assessment) %>% 
#                                 unlist() %>% 
#                                 as.character()
#                       ) %>%
#                       mutate(.pred_Class = factor(.pred_Class))
#                   })
# 
# mod_glm %>% group_by(fold) %>%
#   metrics(truth, .pred_Class)
# # ------------


#------------------------------------------------------------
# create incident likelihood map
#------------------------------------------------------------

# preprocess all input locations
dim(inputs)
baked_all <- bake(recipe, inputs %>% dplyr::select(-train))
dim(baked_all)


# predict incident likelihood at each location
preds_aoi <- predict(mod_final, new_data = baked_all, type = "prob") # get probabilities

preds_aoi <- cbind(preds_aoi, cell_id =baked_all$cell_id)
head(preds_aoi)

# visualize likelihood map
out$likelihood <- NA
out$likelihood[preds_aoi$cell_id] <- preds_aoi$.pred_Incident
plot(out$likelihood, main = 'Charcoaling likelihood')


# incident classification map
preds_aoi <- mutate(preds_aoi, class = ifelse(.pred_Incident < 0.15, 1, 
                                              ifelse(.pred_Incident < 0.5, 2, 
                                                     ifelse(.pred_Incident < 0.85, 3, 4))))
head(preds_aoi)

out$prediction <- NA
out$prediction[preds_aoi$cell_id] <- preds_aoi$class
plot(out$prediction, legend = FALSE, col = c('springgreen4', 'seagreen3', 'indianred1', 'red'), 
     main = 'Charcoaling classification')
legend("bottomleft", legend = c("very unlikely", "unlikely", 'likely', 'very likely'),
       title = 'Charcoaling incident', fill = c('springgreen4', 'seagreen3', 'indianred1', 'red'))

writeRaster(out, paste0('output/model_likelihood_', sc_name, '.tif'), format="GTiff", overwrite = TRUE)


#------------------------------------------------------------
# measure model performance on the test set
#------------------------------------------------------------

test_pred <- preds_aoi %>% filter(cell_id %in% test$cell_id)
test_pred <- cbind(test, dplyr::select(test_pred, -cell_id)) %>%
  dplyr::select(c(cell_id, test, .pred_Unharmed, .pred_Incident, class)) %>% 
  mutate(class = factor(class, levels = c(0,1), labels = c("Unharmed", "Incident")))
head(test_pred)
names(test_pred) <- c('cell_id', 'truth', 'prob_Unharmed', 'prob_Incident', 'prediction')
table(test_pred %>% dplyr::select(- c(cell_id, prob_Unharmed, prob_Incident)), useNA = 'ifany')


# accuracy
metrics(test_pred, truth, prediction)
confusionMatrix(test_pred$prediction, reference = test_pred$truth, positive = 'Incident')
# manually save the information

# ROC curve
roc_df <- dplyr::select(test_pred, c(truth, prob_Incident)) %>% 
  arrange(prob_Incident) %>% 
  mutate(truth = ifelse(truth == 'Incident', 1, 0))
head(roc_df)

# plot ROC curve with AUC value
g <- ggplot(roc_df, aes(m=prob_Incident, d=truth)) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc() +
  geom_abline(intercept = 0, slope = 1, color = 'grey')
g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
ggsave(paste0('output/ROC_curve_', sc_name, '.jpg'))

# alternative: (retrieve tpr and fpr)
# library(pROC)
# ROC <- roc(response = roc_df$truth, predictor = roc_df$prob_Incident, 
#              positive = 'Incident')
# plot(ROC)


#------------------------------------------------------------
# explore features
#------------------------------------------------------------

# dem
par(mfrow = c(2,2))
hist(inputs$elevation_m[inputs$train == 1], breaks = seq(100, 1600, 30), main = 'Incident locations')
hist(inputs$elevation_m[inputs$train == 0], breaks = seq(100, 1600, 30), main = 'Pseudo-absence locations')
hist(inputs$elevation_m, breaks = seq(100, 1600, 30), main = 'Complete area of interest')
par(mfrow = c(1,1))

# dwater
par(mfrow = c(2,2))
hist(inputs$dwater[inputs$train == 1], breaks = seq(0, 50000, 1000), main = 'Incident locations')
hist(inputs$dwater[inputs$train == 0], breaks = seq(0, 50000, 1000), main = 'Pseudo-absence locations')
hist(inputs$dwater, breaks = seq(0, 50000, 1000), main = 'Complete area of interest')
par(mfrow = c(1,1))


hist_f <- function(variable){
  par(mfrow = c(2,2))
  hist(inputs$variable[inputs$train == 1], breaks = seq(100, 1600, 30), main = 'Incident locations')
  hist(inputs$variable[inputs$train == 0], breaks = seq(100, 1600, 30), main = 'Pseudo-absence locations')
  hist(inputs$variable, breaks = seq(100, 1600, 30), main = 'Complete area of interest')
  par(mfrow = c(1,1))
}

hist_f(elevation_m)



# mapping the likelihood and prediction map in leaflet

obscolor <- c("#040404B3","#080918B3","#0E0D24B3","#150F2EB3","#1D1135B3","#24123CB3","#2C1242B3","#341348B3","#3C134EB3","#451353B3","#4D1259B3",
              "#56125DB3","#5F1162B3","#681066B3","#701069B3","#79106DB3","#82106FB3","#8A1172B3","#931373B3","#9B1674B3","#A31A75B3","#AB1E75B3",
              "#B32375B3","#BA2973B3","#C12F71B3","#C8356FB3","#CF3B6BB3","#D64267B3","#DC4962B3","#E2505BB3","#E85752B3","#ED5F48B3","#F2673AB3",
              "#F37133B3","#F47B2CB3","#F58426B3","#F58E23B3","#F69622B3","#F79F25B3","#F7A82CB3","#F7B134B3","#F8B93EB3","#F8C149B3","#F8CA54B3",
              "#F9D25FB3","#F9DB6BB3","#FAE377B3","#FBEC84B3","#FDF490B3","#FFFE9EB3")

pal_lik <- colorNumeric(palette = obscolor, domain = values(out$likelihood),
  na.color = 'transparent')

pal_lik_rev <- colorNumeric(palette = obscolor, domain = values(out$likelihood), 
  reverse = TRUE, na.color = 'transparent')

pal_pred <- colorFactor(palette = c('springgreen4', 'seagreen3', 'indianred1', 'red'),
  domain = values(out$prediction), na.color = 'transparent')

pal_pred(c(1,2,3,4))

library(leaflet)
map <- leaflet() %>% 
  addTiles() %>% 
  addRasterImage(out$likelihood, group = 'likelihood', colors = pal_lik, opacity = 0.75) %>% 
  addLegend(pal = pal_lik_rev, values = values(out$likelihood),
            labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)),
            title = "Probability", group = 'likelihood') %>% 
  addRasterImage(out$prediction, group = 'prediction', opacity = 0.75, 
                 colors = c("#008B45", "#43CD80", "#FF6A6A", "#FF0000")) %>%  # 'springgreen4', 'seagreen3', 'indianred1', 'red'
  addLegend(colors = c("#008B45", "#43CD80", "#FF6A6A", "#FF0000"),
            labels = c('very unlikely', 'unlikely', 'likely', 'very likely'),
            title = "Incident", group = 'prediction') %>%
  addLayersControl(
    overlayGroups = c("likelihood", "prediction"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>% 
  hideGroup('prediction') %>% 
  addScaleBar(
    position = c("bottomleft"),
    options = scaleBarOptions()
  )
map


# JUNK:
# # check relation elevation ~ charcoaling probability
# train$elevation_bin <- ifelse(train$elevation_m>400,
#                           ifelse(train$elevation_m>550,
#                              ifelse(train$elevation_m>700,
#                                   ifelse(train$elevation_m>850,4,3),2),1),0)
# train$incident_num <- 2 - as.numeric(train$train) # check which category got a 0 and which a 1
# 
# aggregate(incident_num~elevation_bin, data=train, FUN=mean)
# 
# preds <- predict(mod_basic, new_data = baked_train, type = "prob")
# preds$truth <- baked_train$train
#   
# plot(baked_train$elevation_m, preds$.pred_Incident)
# plot(baked_train$dwater, preds$.pred_Incident)

