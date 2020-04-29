###### BLOCK 1 ######

library(caret)
library(ggplot2)

# predAndResp_DF is your data with both predictors and response
# list ALL your predictors explicitly in the vector below
# remove them 1 by 1 by running the block and loop below
# and reviewin the results. i'd also copy the top score for each
# run into a seperate document just to review it later
selected_preds2 = c('faith_pd', 'Latitude', 'anti_day_length')

# as you remove them, these lines 
predAndResp_DF2 = predAndResp_DF[,c('Outcome', selected_preds2)]
accuracies = list('full' = log_cm$overall['Accuracy'])

# split the data into two parts 
split <- createDataPartition(y = predAndResp_DF2$Outcome, p = 0.75)
new_train <- predAndResp_DF2[split$Resample1,] 
new_test <- predAndResp_DF2[!(rownames(predAndResp_DF2) %in% rownames(new_train)), ]

# fit the model, predict, and evaluate
log_model <- glm(Outcome~., data=new_train, family = binomial(link="logit"))
log_predict <- predict(log_model, newdata = new_test, type = "response")
log_predict <- ifelse(log_predict > 0.5, levels(predAndResp_DF2$Outcome)[2], levels(predAndResp_DF2$Outcome)[1])
log_cm = confusionMatrix(table(log_predict, new_test$Outcome))

# loop over every predictor 
for (a_col in selected_preds2){
    # drop a col 
    reg_df2short = predAndResp_DF2[,!(colnames(predAndResp_DF2) == a_col)]
    # resplit
    new_train_s <- reg_df2short[split$Resample1,] 
    new_test_s <- reg_df2short[!(rownames(reg_df2short) %in% rownames(new_train_s)), ]
    # redo fit prediction score 
    log_model_s <- glm(Outcome~., data=new_train_s, family = binomial(link="logit"))
    log_predict_s <- predict(log_model_s, newdata = new_test_s, type = "response")
    log_predict_s <- ifelse(log_predict_s > 0.5, levels(c$Outcome)[2], levels(predAndResp_DF2$Outcome)[1])
    # store the results
    accuracies[[a_col]] = confusionMatrix(table(log_predict_s, new_test_s$Outcome))$overall['Accuracy'] - accuracies[['full']]
}

# Print the results
base::sort(unlist(accuracies))

###### BLOCK 2 ######

# this is what you settle on
optimal_model_columns = c('faith_pd', 'anti_day_length', 'Latitude')
# prints a summary of the data in your optimal model
summary(predAndResp_DF[,optimal_model_columns])
# prints the levels of your response
levels(predAndResp_DF$Outcome)

# fit the optimal model on all data
log_model_final <- glm(Outcome~faith_pd + anti_day_length + Latitude, data=predAndResp_DF, 
                       family = binomial(link="logit"))
# print the p-values 
summary(log_model_final)

# this a version of R squared for logistic regressions you can report
# it is the Nagelkerke / Cragg & Uhler version

psuedoR2 <-function(mod, mod0, N) {
    LLf   <- logLik(mod)
    LL0   <- logLik(mod0)
    pR2 = as.vector((1 - exp((2/N) * (LL0 - LLf))) / (1 - exp(LL0)^(2/N)))
    return(pR2)
}

# this is the intercept model
log_model0 <- glm(Outcome~1, data=predAndResp_DF, family = binomial(link="logit"))

psuedoR2(log_model_final, log_model0, nrow(predAndResp_DF))


###### BLOCK 3 ######

# the next two lines make a 300 x 3 matrix containing the median values of each column
# to make a column of a factor at a particular level instead of a median
# use something like factor(rep('Psuedomonas', 100), levels=levels(predAndResp_DF$Genus))

mean_row = apply(predAndResp_DF[,optimal_model_columns], 2, median)
mean_mat = matrix(1, nrow=300,ncol=1) %*% matrix(mean_row, nrow=1, ncol=3)
mean_df = as.data.frame(mean_mat)

# this returns a vector of the 1-100 percentile values of a given column of data
hundred_pts <- function(x) {return(quantile(x, (1:100)/100))}

# this creates a 100 x 3 matrix of each percentile of each column
pct_ranges = apply(predAndResp_DF[,optimal_model_columns], 2, hundred_pts)
colnames(mean_df) <- colnames(pct_ranges)

# this modifies each column to allow each row to vary in one variable  
# while the other two remain at their medians 
mean_df[c(1:100), 'faith_pd'] = pct_ranges[,'faith_pd']
mean_df[c(101:200), 'anti_day_length'] = pct_ranges[,'anti_day_length']
mean_df[c(201:300), 'Latitude'] = pct_ranges[,'Latitude']

# this predicts the response with confidence intervals of the optimal model
model_output <- predict(log_model_final, newdata=mean_df, se.fit = TRUE)[1:2]

# this pulls out a thing to modify the confidence intervals according to a logistic model
# and then dumps the predictions into a data frame
ilink <- family(log_model_final)$linkinv
model_output2 <- data.frame('fit_resp' = ilink(model_output$fit), 
                            'right_upr' = ilink(model_output$fit + (2*model_output$se.fit)),
                            'right_lwr' = ilink(model_output$fit - (2*model_output$se.fit)))

# this creates a single vector of the percentile ranges of each variable in the optimal 
# model 
predictor_range = c(pct_ranges[,'faith_pd'], 
                    pct_ranges[,'anti_day_length'], 
                    pct_ranges[,'Latitude'])

# this creates a vector of the name of the variable in each cell of the previous
var_grps = c(rep('faith_pd', 100), rep('anti_day_length', 100), rep('Latitude', 100))

# this binds those two columns and the predicitons into a single package
mopdf = data.frame('predictor'=predictor_range, 'variable'=var_grps)
mopdf = cbind.data.frame(mopdf, model_output2)
head(mopdf)

# this plots it a line graph with rugs and CI
plt <- ggplot(data = mopdf, aes(x = predictor, y = fit_resp)) + geom_line() +
       geom_rug(aes(x=predictor)) + 
       geom_ribbon(data = mopdf, aes(ymin = right_lwr, ymax = right_upr), alpha = 0.4) + 
       facet_grid(. ~ variable,scales='free')
plt
