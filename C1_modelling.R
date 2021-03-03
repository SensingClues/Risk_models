# Get data splits (train and test set) for the location grid cells.
# Preprocess the data (removing NAs etc.).
# Train the model, also performing cross-validation.
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
# retrieve input locations and corresponding features
#------------------------------------------------------------
inputs <- read.csv('output/location_features.csv')
summary(inputs)
str(inputs)
table(inputs$train)

# correlation check for predictors
corrplot(cor(inputs %>% 
               select(- c(train, test, cell_id, train_label, test_label)) %>% 
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
  fit(train ~ ndvi_max + dwater + droad, data = baked_train)
tidy(mod_final)

# cross-validation
# ----------
set.seed(555)
l_cv <- vfold_cv(baked_train, v = 10, strata = "train")

mod_glm <-map2_df(.x = l_cv$splits,
                  .y = l_cv$id,
                  function (split = .x, fold = .y) 
                  {
                    # Split the data into analysis and assessment tables
                    df_analysis <- analysis(split)
                    df_assessment <- assessment(split)
                      
                    # Build the model
                    mod <-
                      logistic_reg(mode = "classification") %>%
                      set_engine("glm") %>%
                      fit(train ~ ., data = df_analysis)
                      
                    # Summarise Predictions
                    table <- 
                      tibble(fold = fold,
                              truth = df_assessment$train,
                              .pred_Incident = 
                                predict(mod, 
                                        new_data = df_assessment, 
                                        type = "prob")[[".pred_Incident"]],
                              .pred_Unharmed = 
                                predict(mod, 
                                        new_data = df_assessment, 
                                        type = "prob")[[".pred_Unharmed"]],
                              .pred_Class = 
                                predict(mod, 
                                        new_data = df_assessment) %>% 
                                unlist() %>% 
                                as.character()
                      ) %>%
                      mutate(.pred_Class = factor(.pred_Class))
                  })

mod_glm %>% group_by(fold) %>%
  metrics(truth, .pred_Class)
# ------------


#------------------------------------------------------------
# create incident likelihood map
#------------------------------------------------------------
out <- stack('output/raster_template.grd')

# preprocess all input locations
dim(inputs)
baked_all <- bake(recipe, inputs %>% select(-train))
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
preds_aoi <- mutate(preds_aoi, class = ifelse(.pred_Incident > 0.5, 1, 0))

out$prediction <- NA
out$prediction[preds_aoi$cell_id] <- preds_aoi$class
plot(out$prediction, legend = FALSE, col = c('lightblue', 'darkblue'), 
     main = 'Charcoaling classification')
legend("bottomleft", legend = c("not_likely", "likely"),
       title = 'Charcoaling incident', fill = c('lightblue', 'darkblue'))

writeRaster(out, 'output/model_likelihood.tif', format="GTiff", overwrite = TRUE)


#------------------------------------------------------------
# measure model performance on the test set
#------------------------------------------------------------

test_pred <- preds_aoi %>% filter(cell_id %in% test$cell_id)
test_pred <- cbind(test, select(test_pred, -cell_id)) %>%
  select(c(cell_id, test, .pred_Unharmed, .pred_Incident, class)) %>% 
  mutate(class = factor(class, levels = c(0,1), labels = c("Unharmed", "Incident")))
head(test_pred)
names(test_pred) <- c('cell_id', 'truth', 'prob_Unharmed', 'prob_Incident', 'prediction')
table(test_pred %>% select(- c(cell_id, prob_Unharmed, prob_Incident)), useNA = 'ifany')


# accuracy
metrics(test_pred, truth, prediction)
confusionMatrix(test_pred$prediction, reference = test_pred$truth, positive = 'Incident')
# manually save the information

# ROC curve
roc_df <- select(test_pred, c(truth, prob_Incident)) %>% 
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
ggsave('output/ROC_curve.jpg')

# alternative: (retrieve tpr and fpr)
# library(pROC)
# ROC <- roc(response = roc_df$truth, predictor = roc_df$prob_Incident, 
#              positive = 'Incident')
# plot(ROC)


#------------------------------------------------------------
# compare to status quo map
#------------------------------------------------------------

prim <- raster('output/status_quo_likelihood.tif') %>% 
  projectRaster(to = out, method = 'bilinear')
plot(prim)

prim_den <- raster::extract(prim, test$cell_id)
prim_thr <- (range(na.omit(prim_den))[2] - range(na.omit(prim_den))[1]) / 2
prim_pred <- cbind(test, den = prim_den) %>% # check if cbind() works like this
  select(c(cell_id, test, den)) %>% 
  mutate(pred = ifelse(den > prim_thr, 1, 0),
         pred = factor(pred, levels = c(0,1), labels = c("Unharmed", "Incident")))
head(prim_pred)
names(prim_pred) <- c('cell_id', 'truth', 'density', 'prediction')
table(prim_pred %>% select(- c(cell_id, density)), useNA = 'ifany')


# accuracy
metrics(prim_pred, truth, prediction)
confusionMatrix(prim_pred$prediction, reference = prim_pred$truth, positive = 'Incident')
# manually save the information

# number of NAs in status quo predictions on the test set
dim(na.omit(test_pred))[1] - dim(na.omit(prim_pred))[1]


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
