data <- read.csv("project_data/windowed_mean_data_roberta.csv")

data <- data[data$comorbid != "True",]
data <- data[data$gender != "",]
data <- data[data$group_anxiety == 1 | data$group_depression == 1,]

library(lmerTest)
library(tidyverse)
library(splines)
library(ggplot2)

get_hour <- function(date_str){
  minute <- as.numeric(substring(date_str, 15,16))
  hour <- as.numeric(substring(date_str, 12,13))
  return(
    hour + minute / 60
  )
}

data$anxiety <- ifelse(data$group_anxiety == 1,1,0)
data$gender_female <- ifelse(data$gender == "female", 1, 0)
data$weekend <- ifelse(data$day_of_week %in% c(5,6), 1, 0)
data$hour <- data$created_at %>% get_hour
data$hour_sq <- (data$hour - 12) ** 2
data$hour_shifted <- ifelse(data$hour > 20, -(4 + (20 - data$hour)), data$hour)

interval_since_first <- lubridate::interval(lubridate::as_date(min(data$created_at)), data$created_at) %/% lubridate::days(1)
data$days_since_start <- interval_since_first

data <- data %>% arrange(created_at) %>% 
  group_by(user_id) %>% 
  mutate(
    previous_na = lag(ROBERTA_NA),
    time_since_previous_tweet = lubridate::interval(lag(created_at), created_at) %/% lubridate::seconds(1))

data$log_time_since_previous_tweet <- log(data$time_since_previous_tweet)

rmse <- function(y_true, y_hat){
  return(sqrt(mean((y_true - y_hat)^2)))
}

get_model_performance <- function(X, knot_count, degree){
  X <- X[complete.cases(X),]
  print(nrow(X))
  y_true <- c()
  y_hat <- c()
  shuffled_ids <- sample(unique(X$user_id), length(unique(X$user_id)))
  N_splits <- 5
  N_subjects <- length(shuffled_ids)
  N_subjects_per_split <- N_subjects / N_splits
  
  knots <- c()
  knot_step <- 24 / (knot_count + 1)
  for(knot_id in 1:knot_count){
    knots <- c(knots, knot_id * knot_step)
  }
  #knots <- knots - knot_step / 2
  print(knots)
  
  for(cv_index in 0:(N_splits-1)){
    print(cv_index)
    start_index <- 1 + cv_index * floor(N_subjects_per_split)
    end_index <- min(N_subjects, start_index + round(N_subjects_per_split))
    
    sub_ids <- shuffled_ids[start_index:end_index]
    
    #print(start_index)
    #print(end_index)
    #print(length(sub_ids))
    
    join_df <- data.frame(user_id=sub_ids)
    test_df <- X %>% inner_join(join_df, by = "user_id")
    train_df <- X %>% anti_join(join_df, by = "user_id")
    
    #print(nrow(test_df))
    #print(nrow(train_df))
    
    model <- lmer(
      ROBERTA_NA ~ gender_female +
        previous_na +
        log_time_since_previous_tweet +
        anxiety +
        weekend +
        bs(hour,knots = knots,degree = degree) + X19.29 + X30.39 + X..18 + X..40 + (1|user_id)
      ,train_df
    )
    
    y_true <- c(y_true, test_df$ROBERTA_NA)
    y_hat <- c(y_hat, predict(model, test_df, allow.new.levels=T))
  }
  rmse(y_true, y_hat)
}

optimize_model_fit <- function(X, knot_range, degree_range){
  knots <- c()
  degrees <- c()
  rmses <- c()
  
  for(knot_count in knot_range){
    for(degree in degree_range){
      score <- get_model_performance(X, knot_count, degree)
      knots <- c(knots, knot_count)
      degrees <- c(degrees, degree)
      rmses <- c(rmses, score)
    }
  }
  
  data.frame(
    knot_count = knots,
    degree = degrees,
    RMSE = rmses
  )
}

data <- data[is.finite(data$previous_na) & is.finite(data$log_time_since_previous_tweet),]

x <- optimize_model_fit(data, c(2,3,4,5), c(2,3))

model <- lmer(
  ROBERTA_NA ~ gender_female +
    previous_na +
    log_time_since_previous_tweet +
    anxiety +
    weekend +
    bs(hour,knots = c(6, 12, 18),degree = 2) + X19.29 + X30.39 + X..18 + X..40 + (1|user_id)
  ,data
)

summary(model)

### Plot splines

mat <- bs(data$hour,knots = c(6, 12, 18),degree = 2)
head(mat)

fixef(model)[7:11]

Ypred <- mat %*% fixef(model)[7:11] + fixef(model)[1]

Ypred <- scale(Ypred)

ggplot(data.frame(
  hour=data$hour,
  VADER_NA = Ypred
), aes(x=hour, y = VADER_NA))+
  lims(x=c(0,24))+
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  xlab("Time of day (24 hour format)") + 
  ylab("Standardized NA prediction") +
  geom_line()+
  theme_classic()

summary(model)

MuMIn::r.squaredGLMM(model)
MuMIn::r.squaredLR(model)

texreg::texreg(model,)

data_ <- model.frame(model)

data_$residuals <- resid(model)

data_ %>% group_by(anxiety) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )


write.csv(data_, "project_data/residuals_lmm_roberta_6hr.csv")

library(boot)

# Function to compute SD difference
sd_diff <- function(dd, indices) {
  d <- dd[indices, ]
  sd(d$residuals[d$anxiety == 1]) - sd(d$residuals[d$anxiety == 0])
}

# Bootstrap
results <- boot(data_, statistic = sd_diff, R = 1000)
boot.ci(results, type = "perc")

ggplot(data_, aes(group=as.factor(anxiety), fill=as.factor(anxiety), x=residuals))+
  geom_histogram(aes(y=..density..), alpha=0.9)+
  lims(x=c(-1, 1))+
  scale_fill_manual(values=c("#ffbb7c", "#81bcdc"), labels=c("Depression", "Anxiety"), name="Cohort") +
  theme_classic()

car::leveneTest(data_$residuals, as.factor(data_$anxiety), center="median")

data_ %>% group_by(gender_female) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )

ggplot(data_, aes(group=weekend, fill=as.factor(weekend), x=residuals))+
  geom_histogram(aes(y=..density..), alpha=0.9)+
  lims(x=c(-0.4, 0.4))+
  scale_fill_manual(values=c("#ffbb7c", "#81bcdc"), labels=c("Depression", "Anxiety"), name="Cohort") +
  theme_classic()

car::leveneTest(data_$residuals, as.factor(data_$gender_female))

data_ %>% group_by(weekend) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )

car::leveneTest(data_$residuals, as.factor(data_$weekend))

data_ <- data_ %>% mutate(
  previous_tweet_hour_ago = exp(data_$log_time_since_previous_tweet) > 60*60
)

data_ %>% group_by(previous_tweet_hour_ago) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )

car::leveneTest(data_$residuals, as.factor(data_$previous_tweet_hour_ago))

#########
data__ <- data_[!data_$previous_tweet_hour_ago,]

data__ %>% group_by(anxiety) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )

car::leveneTest(data__$residuals, as.factor(data__$anxiety))

data__ %>% group_by(gender_female) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )

car::leveneTest(data__$residuals, as.factor(data__$gender_female))

data__ %>% group_by(weekend) %>% 
  summarize(
    m = mean(residuals),
    sd = sd(residuals)
  )

car::leveneTest(data__$residuals, as.factor(data__$weekend))
