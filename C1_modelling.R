# Get data splits (train and test set) for the location grid cells.
# Preprocess the data (removing NAs etc.).
# Train the model, also performing cross-validation.


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

