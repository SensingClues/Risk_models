# Get data splits (train and test set) for the location grid cells.
# Preprocess the data (removing NAs etc.).
# Train the model, also performing cross-validation.
# Use the final model to create an incident likelihood map.
# Measure model performance on the test set.
# Explore the relationship of the features with the dependent variable (incident yes/no).


# GIS libraries
library(raster)
# tidy libraries
library(tidyverse)
library(tidymodels)
# visualization libraries
library(corrplot)


#------------------------------------------------------------
# retrieve input locations and corresponding features
#------------------------------------------------------------
inputs <- read.csv('output/location_features.csv')
summary(inputs)
str(inputs)
table(inputs$incident)

# correlation check for predictors
corrplot(cor(inputs %>% 
               select(- c(incident, cell_id, cell_label)) %>% 
               filter(! is.na(elevation_m) & ! is.na(slope_deg) & ! is.na(ndvi_mean)))) 


#------------------------------------------------------------
# set divisions
#------------------------------------------------------------

# the training set contains label 0,1 for absence or presence of an indicent
# and label 2 represents the remaining locations in the area of interest
train <- inputs %>% 
  filter(incident != 2) %>% 
  mutate(incident = factor(incident, levels = c(0,1), labels = c("Unharmed", "Charcoal")))
dim(train)
summary(train)
table(train$incident)


# have to retrieve test set still

#------------------------------------------------------------
# preprocessing
#------------------------------------------------------------

recipe <-
  train %>%
  recipe(incident ~ .) %>%
  step_rm(cell_label) %>% 
  step_naomit(all_predictors()) %>% 
  step_corr(all_numeric(), -all_outcomes(), threshold = 0.9) %>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  prep()
# step_normalize(all_numeric(), -all_outcomes()) %>% # grootte coefficenten vergelijken

baked_train <- bake(recipe, train)
dim(baked_train)
summary(baked_train)
table(baked_train$incident)


#------------------------------------------------------------
# training the logistic regression model
#------------------------------------------------------------

mod_basic <- logistic_reg(mode = "classification") %>%
  set_engine("glm") %>%
  fit(incident ~ . -cell_id , data = baked_train) # -dtown -dwater -droad
tidy(mod_basic) 


# step() only accepts glm() model, it gives the same values as the tidymodels model
mod_basic2 <- glm(incident ~ . -cell_id, data = baked_train, family = 'binomial')
summary(mod_basic2)

# stepwise feature selection to determine best predictors
stats::step(mod_basic2) # backwards feature selection


# mod_final <- glm(formula = incident ~ elevation_m + ndvi_mean + ndvi_max, 
#                  family = "binomial", data = baked_train)
# summary(mod_final)
mod_final <- logistic_reg(mode = "classification") %>%
  set_engine("glm") %>%
  fit(incident ~ ndvi_mean + ndvi_max + dwater + droad, data = baked_train)
tidy(mod_final)


# cross-validation
set.seed(555)
l_cv <- vfold_cv(baked_train, v = 10, strata = "incident")

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
                      fit(incident ~ ., data = df_analysis)
                      
                    # Summarise Predictions
                    table <- 
                      tibble(fold = fold,
                              truth = df_assessment$incident,
                              .pred_Charcoal = 
                                predict(mod, 
                                        new_data = df_assessment, 
                                        type = "prob")[[".pred_Charcoal"]],
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


#------------------------------------------------------------
# create incident likelihood map
#------------------------------------------------------------
out <- raster('output/raster_template.grd')

# preprocess all input locations
dim(inputs)
baked_all <- bake(recipe, inputs %>% select(-incident))
dim(baked_all)

# predict charcoaling likelihood at each location
preds_aoi <- predict(mod_final, new_data = baked_all, type = "prob") # get probabilities

preds_aoi <- cbind(preds_aoi, cell_id =baked_all$cell_id)
head(preds_aoi)


# visualize likelihood map
out[] <- NA
out[preds_aoi$cell_id] <- preds_aoi$.pred_Charcoal
plot(out, main = 'Charcoaling likelihood')

# charcoaling classification map
preds_aoi <- mutate(preds_aoi, class = ifelse(.pred_Charcoal > 0.5, 1, 0))

out[preds_aoi$cell_id] <- preds_aoi$class
plot(out, legend = FALSE, col = c('lightblue', 'darkblue'), 
     main = 'Charcoaling classification')
legend("bottomleft", legend = c("not_likely", "likely"),
       title = 'Charcoaling incident', fill = c('lightblue', 'darkblue'))


#------------------------------------------------------------
# measure model performance on the test set
#------------------------------------------------------------

test <- read.csv('output/location_features_test.csv') %>% 
  filter(incident == 1)

test_pred <- cbind(test, pred = preds_aoi$class[test$cell_id]) # check if cbind() works like this
head(test_pred)
names(test_pred) <- c('cell_id', 'truth', 'prediction')
table(test_pred %>% select(-cell_id)) # NAs in pred!
dim(test_pred %>% filter(incident == pred))


#------------------------------------------------------------
# explore features
#------------------------------------------------------------

# dem
par(mfrow = c(2,2))
hist(inputs$elevation_m[inputs$incident == 1], breaks = seq(100, 1600, 30), main = 'Incident locations')
hist(inputs$elevation_m[inputs$incident == 0], breaks = seq(100, 1600, 30), main = 'Pseudo-absence locations')
hist(inputs$elevation_m, breaks = seq(100, 1600, 30), main = 'Complete area of interest')
par(mfrow = c(1,1))

# dwater
par(mfrow = c(2,2))
hist(inputs$dwater[inputs$incident == 1], breaks = seq(0, 50000, 1000), main = 'Incident locations')
hist(inputs$dwater[inputs$incident == 0], breaks = seq(0, 50000, 1000), main = 'Pseudo-absence locations')
hist(inputs$dwater, breaks = seq(0, 50000, 1000), main = 'Complete area of interest')
par(mfrow = c(1,1))


hist_f <- function(variable){
  par(mfrow = c(2,2))
  hist(inputs$variable[inputs$incident == 1], breaks = seq(100, 1600, 30), main = 'Incident locations')
  hist(inputs$variable[inputs$incident == 0], breaks = seq(100, 1600, 30), main = 'Pseudo-absence locations')
  hist(inputs$variable, breaks = seq(100, 1600, 30), main = 'Complete area of interest')
  par(mfrow = c(1,1))
}

hist_f(elevation_m)


# # check relation elevation ~ charcoaling probability
# train$elevation_bin <- ifelse(train$elevation_m>400,
#                           ifelse(train$elevation_m>550,
#                              ifelse(train$elevation_m>700,
#                                   ifelse(train$elevation_m>850,4,3),2),1),0)
# train$incident_num <- 2 - as.numeric(train$incident) # check which category got a 0 and which a 1
# 
# aggregate(incident_num~elevation_bin, data=train, FUN=mean)
# 
# preds <- predict(mod_basic, new_data = baked_train, type = "prob")
# preds$truth <- baked_train$incident
#   
# plot(baked_train$elevation_m, preds$.pred_Charcoal)
# plot(baked_train$dwater, preds$.pred_Charcoal)
