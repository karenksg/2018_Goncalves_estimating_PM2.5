################################################################################
#  Updated version of the R code for the analysis in:
#  GONÇALVES, K.S; WINKLER, MS; BARBOSA-BENCHIMOL, PRB; HOOGH, K; ARTAXO, PA; HACON, SS; CHINDLER, C & KÜNZLI, N. 
#  Development of non-linear models predicting daily fine particle concentrations using aerosol optical depth retrievals and 
#  ground-based measurements at a municipality in the Brazilian Amazon region. 
#  Atmospheric Environment 184 (2018) 156-165. 
#  DOI: https://doi.org/10.1016/j.atmosenv.2018.03.057
#  Last update: June/2021
################################################################################

# LOAD DATA
library(readxl)
modis <- read_excel("C:\\Oberwald\\New Folder\\dataset_validation.txt")

# DEFINING OBJECTS
objects(modis)

aod<-modis$aod
temp<-modis$temp
rh<-modis$rh
aod_t2<-aod*temp^2
aod_t<-aod*temp
aod_rh2<-aod*rh^2
aod_rh<-aod*rh
aod_trh<-aod*temp*rh
aodq_t2<-aod^2*temp^2
aodq_t<-aod^2*temp
aodq_rh2<-aod^2*rh^2
aodq_rh<-aod^2*rh
aodq_trh<-aod^2*temp*rh
aodc_t2<-aod^3*temp^2
aodc_t<-aod^3*temp
aodc_rh2<-aod^3*rh^2
aodc_rh<-aod^3*rh
aodc_trh<-aod^3*temp*rh

rain<-modis$rain

pm25<-modis$pm25
days<-modis$days
cost<-cos(days*2*pi/365.25)
sint<-sin(days*2*pi/365.25)

aodq<-aod^2
aodc<-aod^3
aod_cost<-aod*cost
aod_sint<-aod*sint
aodq_cost<-aodq*cost
aodq_sint<-aodq*sint
aodc_cost<-aodc*cost
aodc_sint<-aodc*sint

####################################################################################################################
# RUNNING LINEAR MODELS
####################################################################################################################
mod1<-lm(pm25~ poly(aod,3) + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh
         + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh)

summary(mod1)
acf(residuals(mod1),type="p")

mod1<-lm(pm25~ poly(aod,3) )
acf(residuals(mod1),type="p")
summary(mod1)



mod1<-lm(pm25~ poly(aod,3) + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh
         + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh + I(aod*rain)+I(aod^2*rain)+I(aod^3*rain))

summary(mod1)
acf(residuals(mod1),type="p")

res<-residuals(mod1)
pred<-predict(mod1)
plot(pred,res)

plot(modis$days,res)

plot(modis$days,pm25)
plot(modis$days,aod)

sum((pm25-pred)^2)

####################################################################################################################
# Try model with log(pm25) and log(aod)
####################################################################################################################

logpm25<-log(pm25)
logaod<-log(aod)

laod_t2<-logaod*temp^2
laod_t<-logaod*temp
laod_rh2<-logaod*rh^2
laod_rh<-logaod*rh
laod_trh<-logaod*temp*rh
laodq_t2<-logaod^2*temp^2
laodq_t<-logaod^2*temp
laodq_rh2<-logaod^2*rh^2
laodq_rh<-logaod^2*rh
laodq_trh<-logaod^2*temp*rh
laodc_t2<-logaod^3*temp^2
laodc_t<-logaod^3*temp
laodc_rh2<-logaod^3*rh^2
laodc_rh<-logaod^3*rh
laodc_trh<-logaod^3*temp*rh


mod2<-lm(pm25~ poly(logaod,3) + laod_t + laod_rh + laod_t2 + laod_rh2 + laod_trh + laodq_t + laodq_rh + laodq_t2 + laodq_rh2 + laodq_trh + laodc_t + laodc_rh + laodc_t2 + laodc_rh2 + laodc_trh + I(logaod*rain) + I(logaod^2*rain))
summary(mod2)

###############################################################################################################
# Comparison of model predictions after excluding period days 300 to 400 with those of from the full model
##############################################################################################################
pm25_<-pm25
pm25_[modis$days>300&modis$days<400]<-NA

mod1b<-lm(pm25_~ poly(aod,3) + aod_t + aod_rh + aod_t2 + aod_rh + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh + I(aod*rain)+I(aod^2*rain)+I(aod^3*rain))
summary(mod1b)
pred_<-predict(mod1b)
pred__<-pred
pred__[modis$days>300 & modis$days<400]<-NA
pred__[!is.na(pred__)]<-pred_
plot(pred__,pred,xlim=c(0,30), ylim=c(0,30))


#################################################################################################################
# Adding sine and cosine terms 
################################################################################################################*


mod1<-lm(pm25~ poly(aod,3) + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh
         + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh + I(aod*rain)+I(aod^2*rain)+I(aod^3*rain)
         + I(aod*cost)+I(aod*sint) + I(aod^2*cost) + I(aod^2*sint) + I(aod^3*cost) + I(aod^3*sint))

summary(mod1)
acf(residuals(mod1),type="p")



mod1<-lm(pm25~ poly(aod,3) + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh 
         + I(aod*cost)+I(aod*sint) + I(aod^2*cost) + I(aod^2*sint) + I(aod^3*cost) + I(aod^3*sint))

pred<-predict(mod1)

mod1b<-lm(pm25_ ~ poly(aod,3) + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + I(aod*cost)+I(aod*sint) + I(aod^2*cost) + I(aod^2*sint) + I(aod^3*cost) + I(aod^3*sint))

pred_<-predict(mod1b)
pred__<-pred
pred__[modis$days>300&modis$days<400]<-NA
pred__[!is.na(pred__)]<-pred_
plot(pred__,pred,xlim=c(0,30), ylim=c(0,30))

#################################################################################################################
# plot of observed and predicted pm25-data
#################################################################################################################

mod1<-lm(pm25~ poly(aod,3))
pred0<-predict(mod1)

plot(days,pm25)
lines(days,pred0,col=c(4))
lines(days,pred,col=c(2))

#################################################################################################################
# Extending the model derived after exclusion of days 300 to 400 to the entire period
#################################################################################################################


mod1b<-lm(pm25_ ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)
summary(mod1b)

newd<-cbind(aod,aodq, aodc,aod_t,aod_rh,aod_t2,aod_rh2,aod_trh,aodq_t,aodq_rh,aodq_t2,aodq_rh2,aodq_trh,
            aod_cost,aod_sint,aodq_cost,aodq_sint,aodq_cost,aodc_sint)

newd<-as.data.frame(newd)

pred<-predict(mod1b,newdata=newd)

plot(days,pred)

#################################################################################################################
### using simplified model
#################################################################################################################

mod1b<-lm(pm25_ ~ aod + aodq + aodc + aod_cost + aod_sint + aodq_cost + aodq_sint)

newd<-cbind(aod,aodq, aodc, aod_cost,aod_sint,aodq_cost,aodq_sint)

newd<-as.data.frame(newd)

pred<-predict(mod1b,newdata=newd)

plot(days,pred)

#######################################################################################################################
# Deriving different predictions after exclusion of values above given thresholds
#######################################################################################################################

mod1_100<-lm(pm25 ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint, subset=(pm25<100))

mod1_90<-lm(pm25 ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint, subset=(pm25<90))

mod1_80<-lm(pm25 ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint, subset=(pm25<80))

mod1_70<-lm(pm25 ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint, subset=(pm25<70))

mod1_60<-lm(pm25 ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint, subset=(pm25<60))

mod1_50<-lm(pm25 ~ aod+aodq+aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aod_cost+ aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint, subset=(pm25<50))

newd<-cbind(aod,aodq, aodc,aod_t,aod_rh,aod_t2,aod_rh2,aod_trh,aodq_t,aodq_rh,aodq_t2,aodq_rh2,aodq_trh,
            aod_cost,aod_sint,aodq_cost,aodq_sint,aodq_cost,aodc_sint)

newd<-as.data.frame(newd)

pred_100<-predict(mod1_100,newdata=newd)
pred_90<-predict(mod1_90,newdata=newd)
pred_80<-predict(mod1_80,newdata=newd)
pred_70<-predict(mod1_70,newdata=newd)
pred_60<-predict(mod1_60,newdata=newd)
pred_50<-predict(mod1_50,newdata=newd)

pairs(~pred_100+pred_90+pred_80+pred_70+pred_60+pred_50)

#######################################################################################################################
# Deriving predictions from model for the log-transformed data
# Problem: backtransformed predictions are too low
# If the residuals were normally distributed, then unbiased predictions would be obtained using the formula
# exp(pred)*exp(s^2/2) where pred are the predictions on the log-scale and s is the standard deviation of residuals on 
# the log-scale.
######################################################################################################################

laod<-log(aod)
laodq<-laod^2
laodc<-laod^3
laod_cost<-laod*cost
laod_sint<-laod*sint
laodq_cost<-laodq*cost
laodq_sint<-laodq*sint
laodc_cost<-laodc*cost
laodc_sint<-laodc*sint

mod4<-lm(logpm25~ laod + laodq+ laodc + laod_t + laod_rh + laod_t2 + laod_rh2 + laod_trh + laodq_t + laodq_rh + laodq_t2 + laodq_rh2 + laodq_trh + laodc_t + laodc_rh + laodc_t2 + laodc_rh2 + laodc_trh  + laod_cost + laod_sint + laodq_cost + laodq_sint + laodc_cost + laodc_sint)

r75=quantile(residuals(mod4),probs=c(0.75))
r25=quantile(residuals(mod4),probs=c(0.25))
sd=(r75-r25)/2/qnorm(0.75)

pred4<-exp(predict(mod4))
pred4corr<-pred4*exp(sd^2/2)

summary(pred4)
summary(pred4corr)
summary(pm25)

# compared with predictions of original model

mod1<-lm(pm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh  + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)

pred1<-predict(mod1)


plot(pred1,pred4corr)
abline(a=0,b=1)

plot(pred1,pred4)
abline(a=0,b=1)

# Sum of squared errors of the different predictions. Clearly show superiority of model 1

sum((pm25-pred1)^2)
sum((pm25-pred4corr)^2)
sum((pm25-pred4)^2)


####### same but after excluding days with pm25 >= 50

mod4<-lm(logpm25~ laod + laodq+ laodc + laod_t + laod_rh + laod_t2 + laod_rh2 + laod_trh + laodq_t + laodq_rh + laodq_t2 + laodq_rh2 + laodq_trh + laodc_t + laodc_rh + laodc_t2 + laodc_rh2 + laodc_trh  + laod_cost + laod_sint + laodq_cost + laodq_sint + laodc_cost + laodc_sint,subset=(pm25<50))

r75=quantile(residuals(mod4),probs=c(0.75))
r25=quantile(residuals(mod4),probs=c(0.25))
sd=(r75-r25)/2/qnorm(0.75)

pred4<-exp(predict(mod4))
pred4corr<-pred4*exp(sd^2/2)

summary(pred4)
summary(pred4corr)
summary(pm25)


mod1<-lm(pm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh  + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint,subset=(pm25<50))

pred1<-predict(mod1)

plot(pred1,pred4corr)
abline(a=0,b=1)

plot(pred1,pred4)
abline(a=0,b=1)

sum((pm25[pm25<50]-pred1)^2)
sum((pm25[pm25<50]-pred4corr)^2)
sum((pm25[pm25<50]-pred4)^2)


#############################################################################################################
# non-linear model:  pm25 = exp(b0 + b1*aod + .,...)
#############################################################################################################

logpm25<-log(pm25)
mod5.<-lm(logpm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)
init<-coefficients(mod5.)
names(init)<-c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13")

mod5<-nls(pm25 ~ exp(b0+b1*aod + b2*aodq+ b3*aodc + b4*aod_t + b5*aod_rh + b6*aod_t2 + b7*aod_rh2 + b8*aod_cost + b9*aod_sint + b10*aodq_cost + b11*aodq_sint + b12*aodc_cost + b13*aodc_sint),
          start=init)
summary(mod5)

pred5<-predict(mod5)
sum((pm25-pred5)^2)

mod1<-lm(pm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh  + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)

pred1<-predict(mod1)
sum((pm25-pred1)^2)

plot(days,pm25)
lines(days,pred1,col=c(4))
lines(days,pred5,col=c(2))

##############################################################################################################
# AR(1) model
##############################################################################################################

mod1<-lm(pm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_trh + aodq_t + aodq_rh + aodq_t2 + aodq_rh2 + aodq_trh + aodc_t + aodc_rh + aodc_t2 + aodc_rh2 + aodc_trh  + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)

summary(mod1)

pred1<-predict(mod1)

sum((pm25-pred1)^2)

pred1<-predict(mod1)
res<-residuals(mod1)
res_1<-c(NA,res[1:length(res)-1])
res_1[1]<-0

mod1b<-update(mod1, ~ . + res_1)
summary(mod1b)

acf(residuals(mod1b),type="p")

pred1b<-predict(mod1b)

sum((pm25-pred1b)^2)

r2<-1-sum((pm25-pred1b)^2)/sum((pm25-mean(pm25))^2)
r2

plot(days,pm25)
lines(days,pred1b,col=c(2))


########### AR(1)-version for non-linear model ################################################################################################


logpm25<-log(pm25)
mod5.<-lm(logpm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)
summary(mod5.)

init<-coefficients(mod5.)
names(init)<-c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13")


mod5b<-nls(pm25 ~ exp(b0+b1*aod + b2*aodq+ b3*aodc + b4*aod_t + b5*aod_rh + b6*aod_t2 + b7*aod_rh2 + b8*aod_cost + b9*aod_sint + b10*aodq_cost + b11*aodq_sint + b12*aodc_cost + b13*aodc_sint ),
           start=init)
summary(mod5b)

res<-residuals(mod5b)
resrel<-res/predict(mod5b)
resrel_1<-c(NA,resrel[1:length(resrel)-1])
resrel_1[1]<-0

init<-c(init,0)

names(init)<-c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14")

mod5b<-nls(pm25 ~ exp(b0+b1*aod + b2*aodq+ b3*aodc + b4*aod_t + b5*aod_rh + b6*aod_t2 + b7*aod_rh2 + b8*aod_cost + b9*aod_sint + b10*aodq_cost + b11*aodq_sint + b12*aodc_cost + b13*aodc_sint+b14*resrel_1), start=init)

acf(residuals(mod5b),type="p")

pred5b<-predict(mod5b)

sum((pm25-pred5b)^2)
r2<-1-sum((pm25-pred5b)^2)/sum((pm25-mean(pm25))^2)
r2 

plot(days,pm25)
lines(days,pred5b,col=c(4))

########## Improvement by adding the square of the lagged residual

logpm25<-log(pm25)
mod5.<-lm(logpm25 ~ aod + aodq+ aodc + aod_t + aod_rh + aod_t2 + aod_rh2 + aod_cost + aod_sint + aodq_cost + aodq_sint + aodc_cost + aodc_sint)

init<-coefficients(mod5.)
names(init)<-c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13")


mod5b<-nls(pm25 ~ exp(b0+b1*aod + b2*aodq+ b3*aodc + b4*aod_t + b5*aod_rh + b6*aod_t2 + b7*aod_rh2 + b8*aod_cost + b9*aod_sint + b10*aodq_cost + b11*aodq_sint + b12*aodc_cost + b13*aodc_sint ),
           start=init)

res<-residuals(mod5b)
resrel<-res/predict(mod5b)
resrel_1<-c(NA,resrel[1:length(resrel)-1])
resrel_1[1]<-0

init<-c(init,0,0)

names(init)<-c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15")

mod5b<-nls(pm25 ~ exp(b0+b1*aod + b2*aodq+ b3*aodc + b4*aod_t + b5*aod_rh + b6*aod_t2 + b7*aod_rh2 + b8*aod_cost + b9*aod_sint + b10*aodq_cost + b11*aodq_sint + b12*aodc_cost + b13*aodc_sint+b14*resrel_1+b15*(resrel_1^2)), start=init)
summary(mod5b)

acf(residuals(mod5b),type="p")

pred5b<-predict(mod5b)

sum((pm25-pred5b)^2)
r2<-1-sum((pm25-pred5b)^2)/sum((pm25-mean(pm25))^2)
r2 

plot(days,pm25)
lines(days,pm25,col=c(2))
lines(days,pred5b,col=c(4))
legend(locator(1),lwd=2.5,lty=1,bty="n", legend=c("PM2.5 measured", "PM2.5 predicted"), col=c("red","blue"))


plot(pm25,pred, ylab="PM2.5 predicted", xlab= "PM2.5 measured",col="dark blue", main="Linear model (Model 1)")
abline(a=0,b=1,col="red")

plot(pm25,pred5b, ylab="PM2.5 predicted", xlab= "PM2.5 measured",col="dark blue", main="Non-linear model (Model 6)")
abline(a=0,b=1,col="red")


####################### Beta predictions ##############

b0<-(2.354e+00)  
b1<-(-1.121e+01)  
b2<-(2.016e+01)  
b3<-(-9.944e+00)  
b4<-(1.706e-01)     
b5<-(5.822e-02)     
b6<-(-3.125e-03)      
b7<-(-4.231e-04)      
b8<-(8.418e+00)  
b9<-(-6.592e+00)  
b10<-(-2.082e+01)  
b11<-(1.009e+01)  
b12<-(9.997e+00)  
b13<-(-4.314e+00)  

PM25<-exp(b0+b1*aod + b2*aodq+ b3*aodc + b4*aod_t + b5*aod_rh + b6*aod_t2 + b7*aod_rh2 + b8*aod_cost + b9*aod_sint + b10*aodq_cost + b11*aodq_sint + b12*aodc_cost + b13*aodc_sint)
PM25<-as.data.frame(PM25)


library (xlsx)
write.xlsx(PM25, "U:/Notes_Backup/pm25_chris.xlsx")