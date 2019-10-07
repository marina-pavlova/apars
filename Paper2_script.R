# This is the paper 2 script file
#
# step 1: restricting analysis to currency crises only
#
cur_data <- read.csv(file="paper_2_currency_only_short.csv", header=TRUE, sep=",")
# checking class bias
table(cur_data$regime_change)
#
# probit models of autocratic regime breakdown during crises
#
# model 1: capital account openness
model1 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_kaopen + kaopen_onset, data=cur_data, family=binomial(link='probit'))
# model 2: monetary independence
model2 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_mi + mi_onset, data=cur_data, family=binomial(link='probit'))
# model 3: exchange rate stability
model3 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers + ers_onset, data=cur_data, family=binomial(link='probit'))
# model 4: capital account openness & monetary independence
model4 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_kaopen + kaopen_onset + d_mi + mi_onset, data=cur_data, family=binomial(link='probit'))
# model 5: capital account openness & exchange rate stability
model5 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_kaopen + kaopen_onset + d_ers + ers_onset, data=cur_data, family=binomial(link='probit'))
# model 6: exchange rate stability & monetary independence
model6 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers + ers_onset + d_mi + mi_onset, data=cur_data, family=binomial(link='probit'))
# model 7: composite index: ers_mi (policy orientation of greater exchange rate stability & monetary independence)
model7 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers_mi + ers_mi_onset, data=cur_data, family=binomial(link='probit'))
# model 8: composite index: mi_kaopen (policy orientation of greater monetary independence and financial openness)
model8 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_mi_kaopen + mi_kaopen_onset, data=cur_data, family=binomial(link='probit'))
# model 9: composite index: ers_kaopen (policy orientation of greater exchange rate stability and financial openness)
model9 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers_kaopen + ers_kaopen_onset, data=cur_data, family=binomial(link='probit'))
#
# creating table 1 (probit)
library(stargazer)
stargazer(model1, model2, model3, model4, model5, model6, model7, model8, model9, type="html", title="Table 1. Adjustment policies and autocratic regime breakdown during crises", out="paper2_table1.doc", dep.var.caption="regime change", dep.var.labels="", keep.stat=c("n", "rsq", "ll"), notes.append=FALSE, notes="Clustered standard errors in brackets. <br> Significance levels: *** p<0.01, ** p<0.05, * p<0.1", notes.align="l")
#
# model diagnostics
#
# checking for multicollinearity
#
library(car)
lapply(list(model1, model2, model3, model4, model5, model6, model7, model8, model9), car::vif)
for (y in 1:9) {
  z=paste0("model",y)
  i=get(z)
  v=paste0("vif",y)
  assign(v, cbind(paste0("model",y), t(car::vif(i))))
}
stargazer(vif1, vif2, vif3, vif4, vif5, vif6, vif7, vif8, vif9, type="html", out="paper2_vif-table.doc", dep.var.caption="VIF", flip=TRUE)
#
# jackknife estimation
#
library(bootstrap)
theta <- function(i, x, dat, coefficient){ coef(glm(i, data = dat[x,], family=binomial(link='probit')))[coefficient] }
for (y in 1:9) {
  z=paste0("model",y)
  i=get(z)
  w=c('d_kaopen', 'kaopen_onset', 'd_mi', 'mi_onset', 'd_ers', 'ers_onset', 'd_kaopen', 'kaopen_onset', 'd_mi', 'mi_onset', 'd_kaopen', 'kaopen_onset', 'd_ers', 'ers_onset', 'd_mi', 'mi_onset', 'd_ers', 'ers_onset', 'd_ers_mi', 'ers_mi_onset', 'd_mi_kaopen', 'mi_kaopen_onset', 'd_ers_kaopen', 'ers_kaopen_onset')
  if (y>=1 & y<=3) { 
    a=getElement(cur_data, w[2*y-1])
    b=getElement(cur_data, w[2*y])
    v=paste0("res",y,"1")
    assign(v, jackknife(1:length(a), theta, i=i, dat=cur_data, coefficient=w[2*y-1]))
    v=paste0("res",y,"2")
    assign(v, jackknife(1:length(b), theta, i=i, dat=cur_data, coefficient=w[2*y]))
  } else if (y>=4 & y<=6) {
    a=getElement(cur_data, w[4*y-9])
    b=getElement(cur_data, w[4*y-8])
    c=getElement(cur_data, w[4*y-7])
    d=getElement(cur_data, w[4*y-6])
    v=paste0("res",y,"1")
    assign(v, jackknife(1:length(a), theta, i=i, dat=cur_data, coefficient=w[4*y-9]))
    v=paste0("res",y,"2")
    assign(v, jackknife(1:length(b), theta, i=i, dat=cur_data, coefficient=w[4*y-8]))
    v=paste0("res",y,"3")
    assign(v, jackknife(1:length(c), theta, i=i, dat=cur_data, coefficient=w[4*y-7]))
    v=paste0("res",y,"4")
    assign(v, jackknife(1:length(d), theta, i=i, dat=cur_data, coefficient=w[4*y-6]))
  } else {
    a=getElement(cur_data, w[2*y+5])
    b=getElement(cur_data, w[2*y+6])
    v=paste0("res",y,"1")
    assign(v, jackknife(1:length(a), theta, i=i, dat=cur_data, coefficient=w[2*y+5]))
    v=paste0("res",y,"2")
    assign(v, jackknife(1:length(b), theta, i=i, dat=cur_data, coefficient=w[2*y+6]))
  }
}
#
# creating data frame with results (only independent variables are included)
#
model_names <- c ('model 1: d_kaopen', 'model 1: kaopen_onset', 'model 2: d_mi', 'model 2: mi_onset', 'model 3: d_ers', 'model 3: ers_onset', 'model 4: d_kaopen', 'model 4: kaopen_onset', 'model 4: d_mi', 'model 4: mi_onset', 'model 5: d_kaopen', 'model 5: kaopen_onset', 'model 5: d_ers', 'model 5: ers_onset', 'model 6: d_mi', 'model 6: mi_onset', 'model 6: d_ers', 'model 6: ers_onset', 'model 7: d_ers_mi', 'model 7: ers_mi_onset', 'model 8: d_mi_kaopen', 'model 8: mi_kaopen_onset', 'model 9: d_ers_kaopen', 'model 9: ers_kaopen_onset')
estimate <- c (as.numeric(model1$coefficients['d_kaopen']), as.numeric(model1$coefficients['kaopen_onset']), as.numeric(model2$coefficients['d_mi']), as.numeric(model2$coefficients['mi_onset']), as.numeric(model3$coefficients['d_ers']), as.numeric(model3$coefficients['ers_onset']), as.numeric(model4$coefficients['d_kaopen']), as.numeric(model4$coefficients['kaopen_onset']), as.numeric(model4$coefficients['d_mi']), as.numeric(model4$coefficients['mi_onset']), as.numeric(model5$coefficients['d_kaopen']), as.numeric(model5$coefficients['kaopen_onset']), as.numeric(model5$coefficients['d_ers']), as.numeric(model5$coefficients['ers_onset']), as.numeric(model6$coefficients['d_mi']), as.numeric(model6$coefficients['mi_onset']), as.numeric(model6$coefficients['d_ers']), as.numeric(model6$coefficients['ers_onset']), as.numeric(model7$coefficients['d_ers_mi']), as.numeric(model7$coefficients['ers_mi_onset']), as.numeric(model8$coefficients['d_mi_kaopen']), as.numeric(model8$coefficients['mi_kaopen_onset']), as.numeric(model9$coefficients['d_ers_kaopen']), as.numeric(model9$coefficients['ers_kaopen_onset']))
st_error <- c (summary(model1)$coefficients['d_kaopen', 2], summary(model1)$coefficients['kaopen_onset', 2], summary(model2)$coefficients['d_mi', 2], summary(model2)$coefficients['mi_onset', 2], summary(model3)$coefficients['d_ers', 2], summary(model3)$coefficients['ers_onset', 2], summary(model4)$coefficients['d_kaopen', 2], summary(model4)$coefficients['kaopen_onset', 2], summary(model4)$coefficients['d_mi', 2], summary(model4)$coefficients['mi_onset', 2], summary(model5)$coefficients['d_kaopen', 2], summary(model5)$coefficients['kaopen_onset', 2], summary(model5)$coefficients['d_ers', 2], summary(model5)$coefficients['ers_onset', 2], summary(model6)$coefficients['d_mi', 2], summary(model6)$coefficients['mi_onset', 2], summary(model6)$coefficients['d_ers', 2], summary(model6)$coefficients['ers_onset', 2], summary(model7)$coefficients['d_ers_mi', 2], summary(model7)$coefficients['ers_mi_onset', 2], summary(model8)$coefficients['d_mi_kaopen', 2], summary(model8)$coefficients['mi_kaopen_onset', 2], summary(model9)$coefficients['d_ers_kaopen', 2], summary(model9)$coefficients['ers_kaopen_onset', 2])
jack_mean <- c (mean(res11$jack.values), mean(res12$jack.values), mean(res21$jack.values), mean(res22$jack.values), mean(res31$jack.values), mean(res32$jack.values), mean(res41$jack.values), mean(res42$jack.values), mean(res43$jack.values), mean(res44$jack.values), mean(res51$jack.values), mean(res52$jack.values), mean(res53$jack.values), mean(res54$jack.values), mean(res61$jack.values), mean(res62$jack.values), mean(res63$jack.values), mean(res64$jack.values), mean(res71$jack.values), mean(res72$jack.values), mean(res81$jack.values), mean(res82$jack.values), mean(res91$jack.values), mean(res92$jack.values))
jack_max <- c (max(res11$jack.values), max(res12$jack.values), max(res21$jack.values), max(res22$jack.values), max(res31$jack.values), max(res32$jack.values), max(res41$jack.values), max(res42$jack.values), max(res43$jack.values), max(res44$jack.values), max(res51$jack.values), max(res52$jack.values), max(res53$jack.values), max(res54$jack.values), max(res61$jack.values), max(res62$jack.values), max(res63$jack.values), max(res64$jack.values), max(res71$jack.values), max(res72$jack.values), max(res81$jack.values), max(res82$jack.values), max(res91$jack.values), max(res92$jack.values))
jack_min <- c (min(res11$jack.values), min(res12$jack.values), min(res21$jack.values), min(res22$jack.values), min(res31$jack.values), min(res32$jack.values), min(res41$jack.values), min(res42$jack.values), min(res43$jack.values), min(res44$jack.values), min(res51$jack.values), min(res52$jack.values), min(res53$jack.values), min(res54$jack.values), min(res61$jack.values), min(res62$jack.values), min(res63$jack.values), min(res64$jack.values), min(res71$jack.values), min(res72$jack.values), min(res81$jack.values), min(res82$jack.values), min(res91$jack.values), min(res92$jack.values))
jack_st_error <- c(res11$jack.se, res12$jack.se, res21$jack.se, res22$jack.se, res31$jack.se, res32$jack.se, res41$jack.se, res42$jack.se, res43$jack.se, res44$jack.se, res51$jack.se, res52$jack.se, res53$jack.se, res54$jack.se, res61$jack.se, res62$jack.se, res63$jack.se, res64$jack.se, res71$jack.se, res72$jack.se, res81$jack.se, res82$jack.se, res91$jack.se, res92$jack.se)
jack_bias <- c(as.numeric(res11$jack.bias), as.numeric(res12$jack.bias), as.numeric(res21$jack.bias), as.numeric(res22$jack.bias), as.numeric(res31$jack.bias), as.numeric(res32$jack.bias), as.numeric(res41$jack.bias), as.numeric(res42$jack.bias), as.numeric(res43$jack.bias), as.numeric(res44$jack.bias), as.numeric(res51$jack.bias), as.numeric(res52$jack.bias), as.numeric(res53$jack.bias), as.numeric(res54$jack.bias), as.numeric(res61$jack.bias), as.numeric(res62$jack.bias), as.numeric(res63$jack.bias), as.numeric(res64$jack.bias), as.numeric(res71$jack.bias), as.numeric(res72$jack.bias), as.numeric(res81$jack.bias), as.numeric(res82$jack.bias), as.numeric(res91$jack.bias), as.numeric(res92$jack.bias))
jack_results <- data.frame(model_names, estimate, st_error, jack_mean, jack_max, jack_min, jack_st_error, jack_bias)
#
# library(knitr)
# kable(jack_results, caption="Jackknife estimation results", format="html)
#
# creating table 7
stargazer(jack_results, title="Table 7. Jackknife estimation results for autocratic regime breakdown during crises", rownames=FALSE, summary=FALSE, type="html", out="paper2_table7.doc")
#
# checking for influential observations
#
library(ggplot2)
outliers_plot <- function(i,y,v) { 
  x=residuals(i, type="partial")
  assign(v, qplot(x, bins=30, main = paste0("Model ", y), xlab="residuals", ylab="count") + theme_bw(), envir = globalenv())
# print(qplot(x, bins=30, main = paste0("Outliers: model ", y), xlab="residuals", ylab="count") + theme_bw())
# Sys.sleep(1)
}
outliers_which <- function(i,x,dat) {
  missing=is.na(x)
  temp.data=subset(dat, !missing)
# print(nrow(temp.data))
  assign("outlier", which(residuals(i, type="partial")==max(residuals(i, type="partial")), arr.ind=TRUE), envir=globalenv())
# print(temp.data[outlier[1],c('cname', 'year')])
# Sys.sleep(1)
}
for (y in 1:9) {
# w=c("d_kaopen", "d_mi", "d_ers", "d_mi_kaopen", "d_ers_kaopen", "d_ers_mi", "d_ers_mi", "d_mi_kaopen", "d_ers_kaopen")
  w=c('d_kaopen', 'kaopen_onset', 'd_mi', 'mi_onset', 'd_ers', 'ers_onset', 'd_kaopen', 'kaopen_onset', 'd_mi', 'mi_onset', 'd_kaopen', 'kaopen_onset', 'd_ers', 'ers_onset', 'd_mi', 'mi_onset', 'd_ers', 'ers_onset', 'd_ers_mi', 'ers_mi_onset', 'd_mi_kaopen', 'mi_kaopen_onset', 'd_ers_kaopen', 'ers_kaopen_onset')
  z=paste0("model",y)
  i=get(z)
  v=paste0("plot",y)
  outliers_plot(i,y,v)
  if (y>=1 & y<=3) {
#   print(z)
    outliers_which(i,getElement(cur_data, w[2*y-1]),cur_data)
    outliers_which(i,getElement(cur_data, w[2*y]),cur_data)
  } else if (y>=4 & y<=6) {
#   print(z)
    outliers_which(i,getElement(cur_data, w[4*y-9]),cur_data)
    outliers_which(i,getElement(cur_data, w[4*y-8]),cur_data)
    outliers_which(i,getElement(cur_data, w[4*y-7]),cur_data)
    outliers_which(i,getElement(cur_data, w[4*y-6]),cur_data)
  } else {
#   print(z)
    outliers_which(i,getElement(cur_data, w[2*y+5]),cur_data)
    outliers_which(i,getElement(cur_data, w[2*y+6]),cur_data)
  }
} 
library(gridExtra)
# gridExtra::grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3)
ggsave("paper2_outliers_plot.png", gridExtra::arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, top="Figure X. Checking for influential outliers"))
#
# removing outliers for models where jackknife bias is large
#
remove_function <- function(y,k,n) {
  v=paste0("res",y,k)
# print(v)
  x=get(v)
  t=get(v)
  m=Mod(as.numeric(x$jack.bias))
  temp.cur_data <- cur_data
  ro=paste0("jack_res",y,k)
  assign(ro, jackknife(1:length(a), theta, i=i, dat=temp.cur_data, coefficient=w[n]), envir=globalenv())
  while (Mod(as.numeric(t$jack.bias))>=m) {
    outliers_which(i,getElement(temp.cur_data,w[n]),temp.cur_data)
    temp.cur_data <- temp.cur_data[-c(outlier[1]),]
    a=getElement(temp.cur_data, w[n])
    vv=paste0("temp")
    assign(vv, jackknife(1:length(a), theta, i=i, dat=temp.cur_data, coefficient=w[n]), envir=globalenv())
    x=get(vv)
    if (Mod(as.numeric(x$jack.bias))<m) {
      t=get(vv)
      ro=paste0("jack_res",y,k)
      assign(ro, jackknife(1:length(a), theta, i=i, dat=temp.cur_data, coefficient=w[n]), envir=globalenv())
    }
    m=Mod(as.numeric(x$jack.bias))
#   print (as.numeric(x$jack.bias))
#   print (as.numeric(t$jack.bias))
  }
}
for (y in 1:9) {
  z=paste0("model",y)
  i=get(z)
  w=c('d_kaopen', 'kaopen_onset', 'd_mi', 'mi_onset', 'd_ers', 'ers_onset', 'd_kaopen', 'kaopen_onset', 'd_mi', 'mi_onset', 'd_kaopen', 'kaopen_onset', 'd_ers', 'ers_onset', 'd_mi', 'mi_onset', 'd_ers', 'ers_onset', 'd_ers_mi', 'ers_mi_onset', 'd_mi_kaopen', 'mi_kaopen_onset', 'd_ers_kaopen', 'ers_kaopen_onset')
    if (y>=1 & y<=3) {
    remove_function(y,1,2*y-1)
    remove_function(y,2,2*y)
  } else if (y>=4 & y<=6) {
    remove_function(y,1,4*y-9)
    remove_function(y,2,4*y-8)
    remove_function(y,3,4*y-7)
    remove_function(y,4,4*y-6)
  } else if (y>=7 & y<=9) {
    remove_function(y,1,2*y+5)
    remove_function(y,2,2*y+6)
  }
}
# creating data frame with new results (after removing outliers)
jack_mean_2 <- c (mean(jack_res11$jack.values), mean(jack_res12$jack.values), mean(jack_res21$jack.values), mean(jack_res22$jack.values), mean(jack_res31$jack.values), mean(jack_res32$jack.values), mean(jack_res41$jack.values), mean(jack_res42$jack.values), mean(jack_res43$jack.values), mean(jack_res44$jack.values), mean(jack_res51$jack.values), mean(jack_res52$jack.values), mean(jack_res53$jack.values), mean(jack_res54$jack.values), mean(jack_res61$jack.values), mean(jack_res62$jack.values), mean(jack_res63$jack.values), mean(jack_res64$jack.values), mean(jack_res71$jack.values), mean(jack_res72$jack.values), mean(jack_res81$jack.values), mean(jack_res82$jack.values), mean(jack_res91$jack.values), mean(jack_res92$jack.values))
jack_max_2 <- c (max(jack_res11$jack.values), max(jack_res12$jack.values), max(jack_res21$jack.values), max(jack_res22$jack.values), max(jack_res31$jack.values), max(jack_res32$jack.values), max(jack_res41$jack.values), max(jack_res42$jack.values), max(jack_res43$jack.values), max(jack_res44$jack.values), max(jack_res51$jack.values), max(jack_res52$jack.values), max(jack_res53$jack.values), max(jack_res54$jack.values), max(jack_res61$jack.values), max(jack_res62$jack.values), max(jack_res63$jack.values), max(jack_res64$jack.values), max(jack_res71$jack.values), max(jack_res72$jack.values), max(jack_res81$jack.values), max(jack_res82$jack.values), max(jack_res91$jack.values), max(jack_res92$jack.values))
jack_min_2 <- c (min(jack_res11$jack.values), min(jack_res12$jack.values), min(jack_res21$jack.values), min(jack_res22$jack.values), min(jack_res31$jack.values), min(jack_res32$jack.values), min(jack_res41$jack.values), min(jack_res42$jack.values), min(jack_res43$jack.values), min(jack_res44$jack.values), min(jack_res51$jack.values), min(jack_res52$jack.values), min(jack_res53$jack.values), min(jack_res54$jack.values), min(jack_res61$jack.values), min(jack_res62$jack.values), min(jack_res63$jack.values), min(jack_res64$jack.values), min(jack_res71$jack.values), min(jack_res72$jack.values), min(jack_res81$jack.values), min(jack_res82$jack.values), min(jack_res91$jack.values), min(jack_res92$jack.values))
jack_st_error_2 <- c(jack_res11$jack.se, jack_res12$jack.se, jack_res21$jack.se, jack_res22$jack.se, jack_res31$jack.se, jack_res32$jack.se, jack_res41$jack.se, jack_res42$jack.se, jack_res43$jack.se, jack_res44$jack.se, jack_res51$jack.se, jack_res52$jack.se, jack_res53$jack.se, jack_res54$jack.se, jack_res61$jack.se, jack_res62$jack.se, jack_res63$jack.se, jack_res64$jack.se, jack_res71$jack.se, jack_res72$jack.se, jack_res81$jack.se, jack_res82$jack.se, jack_res91$jack.se, jack_res92$jack.se)
jack_bias_2 <- c(as.numeric(jack_res11$jack.bias), as.numeric(jack_res12$jack.bias), as.numeric(jack_res21$jack.bias), as.numeric(jack_res22$jack.bias), as.numeric(jack_res31$jack.bias), as.numeric(jack_res32$jack.bias), as.numeric(jack_res41$jack.bias), as.numeric(jack_res42$jack.bias), as.numeric(jack_res43$jack.bias), as.numeric(jack_res44$jack.bias), as.numeric(jack_res51$jack.bias), as.numeric(jack_res52$jack.bias), as.numeric(jack_res53$jack.bias), as.numeric(jack_res54$jack.bias), as.numeric(jack_res61$jack.bias), as.numeric(jack_res62$jack.bias), as.numeric(jack_res63$jack.bias), as.numeric(jack_res64$jack.bias), as.numeric(jack_res71$jack.bias), as.numeric(jack_res72$jack.bias), as.numeric(jack_res81$jack.bias), as.numeric(jack_res82$jack.bias), as.numeric(jack_res91$jack.bias), as.numeric(jack_res92$jack.bias))
jack_results_2 <- data.frame(model_names, estimate, st_error, jack_mean_2, jack_max_2, jack_min_2, jack_st_error_2, jack_bias_2)
#
# creating table X
# stargazer(jack_results_2, title="Table X. Jackknife estimation results after removing outliers", rownames=FALSE, summary=FALSE, type="html", out="paper2_tableX.doc")
#
# bootstrap estimation
library("boot")
bootstrap <- function(formula, data, regressors) {
  dat <- data[regressors,]	
  reg <- glm(formula, data = dat, family=binomial(link='probit')) 
  return(coef(reg)) 
}
for (y in 1:9) {
  z=paste0("model",y)
  w=get(z)
  v=paste0("bs.res",y)
  assign (v, boot(formula=w$formula, data=cur_data, statistic = bootstrap, R=1000))
}
#
# creating data frame with bootstrap results (only independent variables are included)
boot_median <- c (summary(bs.res1)$bootMed[9], summary(bs.res1)$bootMed[10], summary(bs.res2)$bootMed[9], summary(bs.res2)$bootMed[10], summary(bs.res3)$bootMed[9], summary(bs.res3)$bootMed[10], summary(bs.res4)$bootMed[9], summary(bs.res4)$bootMed[10], summary(bs.res4)$bootMed[11], summary(bs.res4)$bootMed[12], summary(bs.res5)$bootMed[9], summary(bs.res5)$bootMed[10], summary(bs.res5)$bootMed[11], summary(bs.res5)$bootMed[12], summary(bs.res6)$bootMed[9], summary(bs.res6)$bootMed[10], summary(bs.res6)$bootMed[11], summary(bs.res6)$bootMed[12], summary(bs.res7)$bootMed[9], summary(bs.res7)$bootMed[10], summary(bs.res8)$bootMed[9], summary(bs.res8)$bootMed[10], summary(bs.res9)$bootMed[9], summary(bs.res9)$bootMed[10])
boot_se <- c (summary(bs.res1)$bootSE[9], summary(bs.res1)$bootSE[10], summary(bs.res2)$bootSE[9],  summary(bs.res2)$bootSE[10], summary(bs.res3)$bootSE[9], summary(bs.res3)$bootSE[10], summary(bs.res4)$bootSE[9], summary(bs.res4)$bootSE[10], summary(bs.res4)$bootSE[11], summary(bs.res4)$bootSE[12], summary(bs.res5)$bootSE[9], summary(bs.res5)$bootSE[10], summary(bs.res5)$bootSE[11], summary(bs.res5)$bootSE[12], summary(bs.res6)$bootSE[9], summary(bs.res6)$bootSE[10], summary(bs.res6)$bootSE[11], summary(bs.res6)$bootSE[12], summary(bs.res7)$bootSE[9], summary(bs.res7)$bootSE[10], summary(bs.res8)$bootSE[9], summary(bs.res8)$bootSE[10], summary(bs.res9)$bootSE[9], summary(bs.res9)$bootSE[10])
boot_bias <- c (summary(bs.res1)$bootBias[9], summary(bs.res1)$bootBias[10], summary(bs.res2)$bootBias[9], summary(bs.res2)$bootBias[10], summary(bs.res3)$bootBias[9], summary(bs.res3)$bootBias[10], summary(bs.res4)$bootBias[9], summary(bs.res4)$bootBias[10], summary(bs.res4)$bootBias[11], summary(bs.res4)$bootBias[12], summary(bs.res5)$bootBias[9], summary(bs.res5)$bootBias[10], summary(bs.res5)$bootBias[11], summary(bs.res5)$bootBias[12], summary(bs.res6)$bootBias[9], summary(bs.res6)$bootBias[10], summary(bs.res6)$bootBias[11], summary(bs.res6)$bootBias[12], summary(bs.res7)$bootBias[9], summary(bs.res7)$bootBias[10], summary(bs.res8)$bootBias[9], summary(bs.res8)$bootBias[10], summary(bs.res9)$bootBias[9], summary(bs.res9)$bootBias[10])
boot_results <- data.frame(model_names, estimate, st_error, boot_median, boot_se, boot_bias)
#
# creating table 8
stargazer(boot_results, title="Table 8. Bootstrap estimation results for autocratic regime breakdown during crises", rownames=FALSE, summary=FALSE, type="html", out="paper2_table8.doc")
#
# running regressions without outliers
library("predictmeans")
for (y in 1:9) {
  z=paste0("model",y)
  w=get(z)
  CookD(w, idn=5)
}
new_cur_data <- cur_data[-c(61, 96, 139),]
model21 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_kaopen + kaopen_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(61, 96),]
model22 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_mi + mi_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(44, 58, 61, 96),]
model23 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers + ers_onset, data=cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(61, 96),]
model24 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_kaopen + kaopen_onset + d_mi + mi_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(44, 61, 96),]
model25 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_kaopen + kaopen_onset + d_ers + ers_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(44, 61, 72, 96),]
model26 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers + ers_onset + d_mi + mi_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(61, 96),]
model27 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers_mi + ers_mi_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(40, 61, 72, 96),]
model28 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_mi_kaopen + mi_kaopen_onset, data=new_cur_data, family=binomial(link='probit'))
new_cur_data <- cur_data[-c(40, 44, 61, 96, 139),]
model29 <- glm(regime_change ~ d_gdppcgr + gdppcgr_onset + vdem_index + duration + party + military + personal + twincrisis + d_ers_kaopen + ers_kaopen_onset, data=new_cur_data, family=binomial(link='probit'))
#
# creating table 6 (probit, without outliers)
stargazer(model21, model22, model23, model24, model25, model26, model27, model28, model29, type="html", title="Table 6. Probit models of autocratic regime breakdown during crises (without outliers)", out="paper2_table6.doc", dep.var.caption="regime change", dep.var.labels="", keep.stat=c("n", "rsq", "ll"), notes.append=FALSE, notes="Clustered standard errors in brackets. <br> Significance levels: *** p<0.01, ** p<0.05, * p<0.1", notes.align="l")
#
#
# step 2: testing hypotheses on the entire population of autocratic country-years
#
full_data <- read.csv(file="paper_2_full_data.csv", header=TRUE, sep=",")
# checking class bias
table(full_data$regime_change)
#
# policy orientation and autocratic regime change
#
# lag variables
#
# library(plyr)
#lg <- function(x)c(NA, x[1:(length(x)-1)])
# full_data2 = ddply(full_data, ~cname, transform, l_vdem=lg(vdem_index))
#
library("stats")
library("pglm")
library("survival")
#
# adjustment policies and autocratic regime change (individual indexes)
#
# model 31: pooled probit - capital account openness & monetary independence
model31 <- glm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + kaopen + kaopen*currency + d_kaopen + d_kaopen*currency + mi + d_mi + d_mi*currency, family=binomial(link="probit"), data=full_data)
# model 32: RE probit - capital account openness & monetary independence
model32 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + kaopen + kaopen*currency + d_kaopen + d_kaopen*currency + mi + d_mi + d_mi*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 33: RE probit with MA independent variables & regional dummies - capital account openness & monetary independence
model33 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + southeastasia + mideast + africa + latam + exussr + europe + currency + kaopen_ma + kaopen_ma*currency + d_kaopen_ma + d_kaopen_ma*currency + mi_ma + d_mi_ma + d_mi_ma*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 34: pooled probit - capital account openness & exchange rate stability
model34 <- glm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + kaopen + kaopen*currency + d_kaopen + d_kaopen*currency + ers + ers*currency + d_ers + d_ers*currency, family=binomial(link="probit"), data=full_data)
# model 35: RE probit - capital account openness & exchange rate stability
model35 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + kaopen + kaopen*currency + d_kaopen + d_kaopen*currency + ers + ers*currency + d_ers + d_ers*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 36: RE probit with MA independent variables & regional dummies - capital account openness & exchange rate stability
model36 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + southeastasia + mideast + africa + latam + exussr + europe + currency + kaopen_ma + kaopen_ma*currency + d_kaopen_ma + d_kaopen_ma*currency + ers_ma + ers_ma*currency + d_ers_ma + d_ers_ma*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 37: pooled probit - exchange rate stability & monetary independence
model37 <- glm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + mi + d_mi + d_mi*currency + ers + ers*currency + d_ers + d_ers*currency, family=binomial(link="probit"), data=full_data)
# model 38: RE probit - exchange rate stability & monetary independence
model38 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + mi + d_mi + d_mi*currency + ers + ers*currency + d_ers + d_ers*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 39: RE probit with MA independent variables & regional dummies - exchange rate stability & monetary independence
model39 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + southeastasia + mideast + africa + latam + exussr + europe + currency + mi_ma + d_mi_ma + d_mi_ma*currency + ers_ma + ers_ma*currency + d_ers_ma + d_ers_ma*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
#
# policy orientation and autocratic regime change (composite indexes)
# 
# model 41: pooled probit - ers_mi (policy orientation of greater exchange rate stability and monetary independence)
model41 <- glm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + ers_mi + d_ers_mi + d_ers_mi*currency, family=binomial(link="probit"), data=full_data)
# model 42: RE probit - ers_mi (policy orientation of greater exchange rate stability and monetary independence)
model42 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + ers_mi + d_ers_mi + d_ers_mi*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 43: RE probit with MA independent variables & regional dummies - ers_mi (policy orientation of greater exchange rate stability and monetary independence)
model43 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + southeastasia + mideast + africa + latam + exussr + europe + currency + ers_mi_ma + d_ers_mi_ma + d_ers_mi_ma*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 44: pooled probit - mi_kaopen (policy orientation of greater monetary independence and financial openness)
model44 <- glm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + mi_kaopen + d_mi_kaopen + d_mi_kaopen*currency, family=binomial(link="probit"), data=full_data)
# model 45: RE probit - mi_kaopen (policy orientation of greater monetary independence and financial openness)
model45 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + mi_kaopen + d_mi_kaopen + d_mi_kaopen*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 46: RE probit with MA independent variables & regional dummies - mi_kaopen (policy orientation of greater monetary independence and financial openness)
model46 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + southeastasia + mideast + africa + latam + exussr + europe + currency + mi_kaopen_ma + d_mi_kaopen_ma + d_mi_kaopen_ma*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 47: pooled probit - ers_kaopen (policy orientation of greater exchange rate stability and financial openness)
model47 <- glm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + ers_kaopen + d_ers_kaopen + d_ers_kaopen*currency, family=binomial(link="probit"), data=full_data)
# model 48: RE probit - ers_kaopen (policy orientation of greater exchange rate stability and financial openness)
model48 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + currency + ers_kaopen + d_ers_kaopen + d_ers_kaopen*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
# model 49: RE probit with MA independent variables & regional dummies - ers_kaopen (policy orientation of greater exchange rate stability and financial openness)
model49 <- pglm (regime_change ~ gdppcgr + vdem_index + gwf_duration + prevrc + wr_party + wr_military + wr_personal + southeastasia + mideast + africa + latam + exussr + europe + currency + ers_kaopen_ma + d_ers_kaopen_ma + d_ers_kaopen_ma*currency, model="random", family=binomial("probit"), index=c("cname", "year"), data=full_data)
#
# creating tables 2 and 3
#
library(texreg)
extract.pglm <- function (model, include.nobs = TRUE, include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$estimate)
  coefficients <- s$estimate[, 1]
  standard.errors <- s$estimate[, 2]
  significance <- s$estimate[, 4]
  loglik.value <- s$loglik
  n <- nrow(model$model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, loglik.value)
    gof.names <- c(gof.names, "Log-Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                     se = standard.errors, pvalues = significance, gof.names = gof.names, 
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)
}
setMethod("extract", signature = className("maxLik", "maxLik"), definition = extract.pglm)
#
# creating table 2
library(texreg)
htmlreg(list(model31, model32, model33, model34, model35, model36, model37, model38, model39), file="paper2_table2.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), custom.coef.names=c(NA, NA, NA, "duration", NA, "party", "military", "personal", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "kaopen", "d_kaopen", "mi", "d_mi", "currency:kaopen", "currency:d_kaopen", "currency:d_mi", NA, NA, NA, NA, "ers", "d_ers", "currency:ers", "currency:d_ers"), reorder.coef=c(2, 3, 4, 5, 6, 7, 8, 18, 19, 20, 21, 22, 23, 17, 1, 9, 10, 14, 11, 15, 12, 13, 16, 24, 26, 25, 27), caption.above=TRUE, caption="<b>Table 2. Adjustment policies and autocratic regime breakdown</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: *** p<0.01, ** p<0.05, * p<0.1", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE)
#
# creating table 3
htmlreg(list(model41, model42, model43, model44, model45, model46, model47, model48, model49), file="paper2_table3.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), custom.coef.names=c(NA, NA, NA, "duration", NA, "party", "military", "personal", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "ers_mi", "d_ers_mi", "currency:d_ers_mi", NA, NA, NA, "mi_kaopen", "d_mi_kaopen", "currency:d_mi_kaopen", NA, NA, NA, "ers_kaopen", "d_ers_kaopen", "currency:d_ers_kaopen"), reorder.coef=c(2, 3, 4, 5, 6, 7, 8, 14, 15, 16, 17, 18, 19, 13, 1, 9, 10, 11, 12, 20, 21, 22, 23, 24, 25), caption.above=TRUE, caption="<b>Table 3. Policy orientation and autocratic regime breakdown</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: *** p<0.01, ** p<0.05, * p<0.1", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE)
#
# model diagnostics
#
# checking for multicollinearity
library(car)
lapply(list(model31, model34, model37, model41, model44, model47), car::vif)
#
# creating figures
#
library("jtools")
effect_plot(model21, pred=d_kaopen, interval = TRUE, plot.points = TRUE, data=cur_data)
plot_summs(model1, model2, model3, scale = TRUE)
#
# appendix
#
# descriptive statistics
#
library(psych)
des1 = describe(cur_data, fast=TRUE, omit=TRUE)
des1 <- des1[-c(1, 2, 3, 4, 5, 6, 9, 10, 11, 13, 18, 20, 21, 25, 26, 27, 28, 29, 30, 33, 34, 48, 49, 50, 58, 59), -c(1, 7, 8)]
print(des1, digits=3)
des2 = describe(full_data, fast=TRUE, omit=TRUE)
des2 <- des2[-c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 44, 46, 47, 48, 50, 51, 57, 74), -c(1, 7, 8)]
print(des2, digits=3)
library(stargazer)
# creating tables 4 and 5
stargazer(des1, title="Table 4. Descriptive statistics for 'crisis only' data", summary=FALSE, type="html", out="paper2_table4.doc")
stargazer(des2, title="Table 5. Descriptive statistics for panel data", summary=FALSE, type="html", out="paper2_table5.doc")
