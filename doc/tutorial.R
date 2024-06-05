## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5, fig.height=3
  ) 
library(knitr)

## -----------------------------------------------------------------------------
  library(predtools)
  data(gusto)
  set.seed(1)

## -----------------------------------------------------------------------------
  gusto$y <- gusto$day30
  gusto$kill <- (as.numeric(gusto$Killip)>1)*1

## -----------------------------------------------------------------------------
  dev_data <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),] #The regl variable contains location codes
  val_data <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]
  model <- glm(y ~ age + miloc + pmi + kill + pmin(sysbp,100) + pulse, data=dev_data, family=binomial(link="logit"))

## ----echo=FALSE---------------------------------------------------------------
  kable(cbind("Coefficients"=summary(model)$coefficients[,1]))

## -----------------------------------------------------------------------------
  val_data$pi  <- predict(model, type="response", newdata=val_data)

## ----fig.width = 4, fig.height = 4--------------------------------------------
  predtools::calibration_plot(val_data, obs="y", pred="pi")

## -----------------------------------------------------------------------------
  library(cumulcalib)
  res <- cumulcalib(val_data$y, val_data$pi) 

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
  plot(res)

## -----------------------------------------------------------------------------
  val_data$pi2 <- val_data$pi*0.75/(1-val_data$pi*(1-0.75)) #One-shot transformation of risk to odds and back
  res <- cumulcalib(val_data$y, val_data$pi2) 
  predtools::calibration_plot(val_data, obs="y", pred="pi2")
  plot(res)

## -----------------------------------------------------------------------------
  val_data$pi2 <- val_data$pi*1.25/(1-val_data$pi*(1-1.25)) #One-shot transformation of risk to odds and back
  res <- cumulcalib(val_data$y, val_data$pi2) 
  predtools::calibration_plot(val_data, obs="y", pred="pi2")
  plot(res)

## -----------------------------------------------------------------------------
  dev_data2 <- dev_data[sample(nrow(dev_data), 500, replace=F),]
  model2 <- glm(y ~ age + miloc + pmi + kill + pmin(sysbp,100) + pulse, data=dev_data2, family=binomial(link="logit"))
  val_data$pi2  <- predict(model2, type="response", newdata=val_data)

## -----------------------------------------------------------------------------
  predtools::calibration_plot(val_data, obs="y", pred="pi2")

## -----------------------------------------------------------------------------
  res2 <- cumulcalib(val_data$y, val_data$pi2)
  summary(res2)
  plot(res2)

