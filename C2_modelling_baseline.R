# Measure performance of the status quo likelihood map (script A3)
# on the independent test set. 
# For comparison with the newly created model (script C1), the area
# that is NA in the status quo likelihood map is filled in with
# values randomnly sampled from the existing density values, allowing
# to create a ROC curve.


# GIS libraries
library(raster)
library(rgdal)
# modelling libraries
library(tidyverse)
library(caret)
# visualization libraries
library(plotROC)


#------------------------------------------------------------
# setup
#------------------------------------------------------------

# choose a scenario name
sc_name <- 'SC_1month' # other scenarios: 'SC_6month', 'SC_3month', 'SC_1month'

inputs <- read.csv(paste0('output/location_features_', sc_name, '.csv'))
test <- inputs %>% 
  filter(test != 2) %>% 
  mutate(test = factor(test, levels = c(0,1), labels = c("Unharmed", "Incident")))
train <- inputs %>% 
  filter(train != 2) %>% 
  mutate(train = factor(train, levels = c(0,1), labels = c("Unharmed", "Incident")))

out <- stack(paste0('output/rasters_', sc_name, '.tif'))


#------------------------------------------------------------
# baseline: random assignment corrected for incident probability
#------------------------------------------------------------

# using rnorm assumes normal distribution ???
# random_lik <- rnorm(length(inputs$cell_id), mean = 0.5, sd = 0.15)
# 
# out$base_lik <- NA
# out$base_lik[inputs$cell_id] <- random_lik
# plot(out$base_lik)

# incident classification map
# random_class <- ifelse(random_lik > 0.5, 1, 0)
incident_prob <- length(filter(train, train_label == 'incident')$train)/length(inputs$cell_id)
set.seed(25)
random_class <- sample(c(0,1), size = length(inputs$cell_id), replace = TRUE, 
                       prob = c(1-incident_prob, incident_prob))

out$base_pred <- NA
out$base_pred[inputs$cell_id] <- random_class
plot(out$base_pred, legend = FALSE, col = c('lightblue', 'darkblue'), 
     main = 'Charcoaling classification')
legend("bottomleft", legend = c("not_likely", "likely"),
       title = 'Charcoaling incident', fill = c('lightblue', 'darkblue'))

base_pred <- raster::extract(out$base_pred, test$cell_id)
test_pred <- cbind(test, base_pred) %>%
  dplyr::select(c(cell_id, test, base_pred)) %>% 
  mutate(base_pred = factor(base_pred, levels = c(0,1), labels = c("Unharmed", "Incident")))
names(test_pred) <- c('cell_id', 'truth', 'prediction')
table(test_pred %>% dplyr::select(-cell_id), useNA = 'ifany')

confusionMatrix(test_pred$prediction, reference = test_pred$truth, positive = 'Incident')

# plot ROC curve with AUC value
roc_df <- as.data.frame(rep(0.5, length(test$cell_id)))
roc_df <- cbind(roc_df, test$test)
names(roc_df) <- c('prob_Incident', 'truth')
roc_df <- mutate(roc_df, truth = ifelse(truth == 'Incident', 1, 0))

g <- ggplot(roc_df, aes(m=prob_Incident, d=truth)) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc() +
  geom_abline(intercept = 0, slope = 1, color = 'grey')
g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", 0.5))


#------------------------------------------------------------
# status quo likelihood: incident density corrected for patrolling density - with NA
#------------------------------------------------------------

prim <- raster(paste0('output/status_quo_likelihood_', sc_name,'.tif')) %>% 
  projectRaster(to = out, method = 'bilinear')
plot(prim)

prim_den <- raster::extract(prim, test$cell_id)
prim_thr <- (range(na.omit(prim_den))[2] - range(na.omit(prim_den))[1]) / 2
prim_pred <- cbind(test, den = prim_den) %>% # check if cbind() works like this
  dplyr::select(c(cell_id, test, den)) %>% 
  mutate(pred = ifelse(den > prim_thr, 1, 0),
         pred = factor(pred, levels = c(0,1), labels = c("Unharmed", "Incident")))
names(prim_pred) <- c('cell_id', 'truth', 'density', 'prediction')
head(prim_pred)
table(prim_pred %>% dplyr::select(- c(cell_id, density)), useNA = 'ifany')


# accuracy
metrics(prim_pred, truth, prediction)
confusionMatrix(prim_pred$prediction, reference = prim_pred$truth, positive = 'Incident')
# manually save the information

# number of NAs in status quo predictions on the test set
dim(na.omit(test_pred))[1] - dim(na.omit(prim_pred))[1]

# ROC not possible because of large amount of NA


#------------------------------------------------------------
# status quo likelihood: incident density corrected for patrolling density - without NA
#------------------------------------------------------------

# retrieve the area of interest polygon
aoi <- readOGR(dsn = "data", layer = "study_area") %>% 
  spTransform(crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs"))

# create a buffer around the initial status quo likelihood map
# (because the edges are abrupt but at borders the densities are probably similar)
prim_buf <- focal(prim, w=matrix(1,3,3), fun=mean, NAonly=TRUE, na.rm=TRUE)
plot(prim_buf)
plot(aoi, add=TRUE)


# for filling in the rest of the AOI:

# get the initial density values to get distribution parameters mean and sd
dens_vals <- na.omit(values(prim))
hist(dens_vals, breaks = seq(0, 12.5, 0.5))

# # to take a look at the distribution of the status quo density values
# dens_vals <- sort(dens_vals)
# x <- seq(-4,4,length=1000)*sqrt(var(dens_vals)) + mean(dens_vals)
# y <- dnorm(x,mean=mean(dens_vals),sd=sqrt(var(dens_vals)))
# plot(y~x, axes=TRUE, ylab="",xlab="",type="l", col="blue",lwd=3) # ,xlim=xlim,ylim=range(y)/V1size

# get the cell numbers that have values after the buffer was applied
prim_buf_df <- cellsFromExtent(prim_buf, extent(prim_buf)) %>% 
  cbind(values(prim_buf)) %>% 
  as.data.frame() %>% 
  dplyr::filter(! is.na(V2))
names(prim_buf_df) <- c('cell_id', 'vals')

# get the cell numbers of the AOI that are NA
fill_cell_ids <- setdiff(inputs$cell_id, prim_buf_df$cell_id)

# draw density values from the distribution of initial density values
set.seed(25)
fill_vals <- sample(dens_vals, size = length(fill_cell_ids), replace = TRUE)
# fill_vals <- rnorm(n = length(fill_cell_ids), # length(inputs$cell_id) - length(na.omit(values(prim_buf)))
#                    mean=mean(dens_vals),sd=sqrt(var(dens_vals)))
# # density cannot be below zero, artifact of normal distribution, set values below zero to zero
# fill_vals <- ifelse(fill_vals < 0, 0, fill_vals)
hist(fill_vals, breaks = seq(0, 12.5, 0.5))


# add the AOI density values to the status quo with buffer raster
prim_buf[fill_cell_ids] <- fill_vals
plot(prim_buf)
plot(mask(prim_buf, aoi))

# turn density values into probabilities
prim_buf_lik <- prim_buf / max(na.omit(values(prim_buf)))
plot(prim_buf_lik)
plot(mask(prim_buf_lik,aoi), main = 'Status quo likelihood')

# get yes/no incident prediction based on likelihood
prim_buf_prob <- raster::extract(prim_buf_lik, test$cell_id)

# predict using probabilities
prim_buf_pred <- cbind(test, prob = prim_buf_prob) %>% # check if cbind() works like this
  dplyr::select(c(cell_id, test, prob)) %>% 
  mutate(pred = ifelse(prob > 0.5, 1, 0),
         pred = factor(pred, levels = c(0,1), labels = c("Unharmed", "Incident")))
names(prim_buf_pred) <- c('cell_id', 'truth','prob_Incident', 'prediction')
head(prim_buf_pred)
table(prim_buf_pred %>% dplyr::select(- c(cell_id, prob_Incident)), useNA = 'ifany')


# # predict using densities and a threshold at density range midpoint
# prim_buf_den <- raster::extract(prim_buf, test$cell_id)
# prim_buf_thr <- (range(na.omit(values(prim_buf)))[2] - range(na.omit(values(prim_buf)))[1]) / 2
# prim_buf_pred <- cbind(test, den = prim_buf_den) %>% # check if cbind() works like this
#   dplyr::select(c(cell_id, test, den)) %>% 
#   mutate(pred = ifelse(den > prim_buf_thr, 1, 0),
#          pred = factor(pred, levels = c(0,1), labels = c("Unharmed", "Incident")))
# names(prim_buf_pred) <- c('cell_id', 'truth', 'density', 'prediction')
# head(prim_buf_pred)
# table(prim_buf_pred %>% dplyr::select(- c(cell_id, density)), useNA = 'ifany')


# accuracy
metrics(prim_buf_pred, truth, prediction)
confusionMatrix(prim_buf_pred$prediction, reference = prim_buf_pred$truth, positive = 'Incident')
# manually save the information

# number of NAs in status quo predictions on the test set
dim(na.omit(test_pred))[1] - dim(na.omit(prim_pred))[1]


# ROC curve
roc_df <- dplyr::select(prim_buf_pred, c(truth, prob_Incident)) %>% 
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
ggsave(paste0('output/status_quo_ROC_curve_', sc_name, '.jpg'))

