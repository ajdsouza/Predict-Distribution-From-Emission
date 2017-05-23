#
# author : ajdsouza31 
#
#  Implement data analysis and model fitting a fit a model
# to predict the mean and variance of a random variable that
# depends on a function of other random variables which
# have different distributions
#
#####

set.seed(20160227) ### set the random seed


library(lattice)
library(glmnet)
library(corrplot)
library(GGally)

#---------------------------------------
# Functions
#---------------------------------------


#-------------------------------------------------------------------
# X1 vs X2 with display data
#-------------------------------------------------------------------
scatterplot.xy <- function(file.name,
	data.xy,
	display.data,
	title.label,
	title.main ) {

	png(paste(file.name,"png",sep="."),width=1000,height=800)

	m <- matrix(c(1,2,2),nrow = 3,ncol = 1,byrow = TRUE)

	par(oma=c(1,1,1,1))

	layout(mat = m,heights = c(0.80,0.1,.1))

	par(mar = c(4,4,1,0))

	plot(data.xy$X1,data.xy$X2,type='p', pch=21, col=as.numeric(display.data), 
		bg=as.numeric(display.data),
		xlab="X1", ylab="X2")

	title(main=title.main, col.main='red', font.main=4,outer=FALSE)

	par(mar = c(0,0,0,0))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="")

	legend(x="bottom", inset=0, paste(title.label,levels(as.factor(Vqt))),
		 fill=levels(as.factor(display.data)) ,cex=1, horiz = TRUE)

	dev.off()

}


#-----------------------------------------------------------
# X vs Residuals, by display data
#----------------------------------------------------------
plotxy.residuals <- function (
	file.name,
	data,
	model.residuals,
	std.residuals,
	display.data,
	true.value,
	fitted.value,
	title.fitted,
	title.legend,
	title.model
	) {

	png(paste(file.name,"png",sep="."),width=1000,height=800)

	m <- matrix(c(1,2,3,4,5,5),nrow = 6,ncol = 1,byrow = TRUE)

	par(oma=c(1,1,1,1))

	layout(mat = m,heights = c(0.050,0.3,0.3,0.3,0.025,0.025))

	# title
	par(mar = c(0,0,0,4))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="")
	##title(main=paste("Fitted/Residuals ",title.model,title.fitted,sep=" - "), 
	##	col.main='red', font.main=4, cex=1.5, outer=FALSE)


	# fitted vs true values
	par(mar = c(4,4,1,0))
	plot(true.value,fitted.value,
		col=as.numeric(display.data),bg=as.numeric(display.data),
		type='p', pch=21,
		xlab=title.fitted, ylab=paste("Fitted - ",title.fitted," - ",title.model,sep="") )

	title(main=paste("Fitted Vs True Value",title.fitted,sep=""),
		col.main='blue', font.main=4,outer=FALSE)


	# a scatter plot of Mean(Y) vs Standardized Residuals
	par(mar = c(4,4,1,0))
	plot(data$muhat,std.residuals,
		col=as.numeric(display.data),bg=as.numeric(display.data),
		type='p', pch=21,
		xlab="Mean(Y)", 
		ylab=paste("Standardized Residuals for ",title.fitted,sep=""))
	
	title(main=paste("Mean(Y) Vs Standardized Residuals by ",title.legend,sep=""),
		col.main='blue', font.main=4,outer=FALSE)


	# a scatter plot of Variance(Y) vs Standardized Residuals
	par(mar = c(4,4,1,0))
	plot(data$Vhat,std.residuals,
		col=as.numeric(display.data),bg=as.numeric(display.data),
		type='p', pch=21,
		xlab="Variance(Y)", 
		ylab=paste("Standardized Residuals for ",title.fitted,sep=""))
	
	title(main=paste("Variance(Y) Vs Standardized Residuals by ",title.legend,sep=""),
		col.main='blue', font.main=4,outer=FALSE)


	# legend
	par(mar = c(0,0,0,0))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="")
	legend(x="bottom", inset=0, paste(title.legend,levels(as.factor(muqt))),
		 fill=levels(as.factor(muqt)) ,cex=1, horiz = TRUE)


	dev.off()

}














#--------------------------------------------------------------
#  Read the data
#---------------------------------------------------------------

## Read Training Data
distpredtrain <- read.table(file = "http://www2.isye.gatech.edu/~ymei/7406/midtermtrain.csv", sep=",")


## Testing Data
distpredtest  <- read.table(file = "http://www2.isye.gatech.edu/~ymei/7406/midtermtest.csv", sep=",")


#----------------------------------------------------------------------
# Data preparation
#----------------------------------------------------------------------
## Some plots for exploratory data analysis
X1 <- distpredtrain[,1]
X2 <- distpredtrain[,2]
muhat <- apply(distpredtrain[,3:202], 1, mean)
Vhat  <- apply(distpredtrain[,3:202], 1, var)


## regression with poly terms
poly.x1.max <- 6
poly.x2.max <- 6

X1_poly <- poly(X1,poly.x1.max,raw=TRUE)[,-1]
X2_poly <- poly(X2,poly.x2.max,raw=TRUE)[,-1]

data0 <- data.frame(X1 = X1, 
	X2=X2,
	X1_poly,
	X2_poly,
	muhat = muhat, 
	Vhat = Vhat)

# muhat and vhat quartiles
muqt <- as.integer(cut(muhat, quantile(muhat, probs=0:4/4), include.lowest=TRUE))
Vqt <- as.integer(cut(Vhat, quantile(Vhat, probs=0:4/4), include.lowest=TRUE))



#---------------------------------------------------------------------
# Data Exploration and Analysis
#---------------------------------------------------------------------
## dim=2911x202 
## The first two columns are X1 and X2 values, and the last 200 columns are the Y valus
dim(distpredtrain)

## This should be a 1066*2 matrix
## Please add two columns for your estimation of the mean and variance of the Y variable. 
dim(distpredtest)


#---------------------------------------------------------------------
# Box and Density Plots - Data Exploration and Analysis
#----------------------------------------------------------------------


png("mt_train_data_anal_1.png",width=1000,height=800)

m <- matrix(c(1,2,3,4),nrow = 4,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(0.05,0.35,0.3,0.3))

par(mar = c(0,0,0,5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main='Training Data Analysis', col.main='red', font.main=4, cex=1.5, outer=TRUE)

par(mar = c(4,4,1,0))

boxplot(cbind(X1=distpredtrain[,1], X2=distpredtrain[,2],Y=stack(distpredtrain[,3:202])[,1],
		"Mean(Y)"=muhat,"Variance(Y)"=Vhat),
		col=(c("gold","darkgreen")),
		main="Boxplot X1,X2,Y, Mean(Y) and Variance(Y) ",
		col.main='red', font.main=4)

par(mar = c(4,4,1,0))

plot(density(muhat), main="Density Plot - Main(Y)")
polygon(density(muhat),col="red", border="blue")

par(mar = c(4,4,1,0))

plot(density(Vhat), main="Density Plot - Variance(Y)")
polygon(density(Vhat),col="red", border="blue")

dev.off()






png("mt_train_data_anal_density_plot.png",width=1000,height=800)

m <- matrix(c(1,2,3,4),nrow = 4,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(0.25,0.25,0.25,0.25))


for (i in c(100,1000,1500,2500)) {
	par(mar = c(4,4,1,0))

	y.density <- density(stack(distpredtrain[i,3:202])[,1])

	plot(y.density, main=paste("Density Plot of Y for X1=",distpredtrain[i,1],"X2=",distpredtrain[i,2]))

	polygon(y.density,col="red", border="blue")
}

dev.off()

#---------------------------------------------------------------------
# Correlation Plots - Data Exploration and Analysis
#----------------------------------------------------------------------
# The corrlation table (the last column is Y)
corr=round(cor( data0[,c("X1","X2","muhat","Vhat")] ),2)
corr
png("mt_corr_plot.png",width=1000,height=800)

m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(0.5,0.5))

par(mar = c(4,4,1,0))

corrplot(corr, order = "AOE", cl.ratio = 0.2, cl.align = "r",
		tl.pos = "d",tl.srt = 60)
title(main="Correlation Plot - X1,X2,Mean,Variance",
	col.main='red', font.main=4,outer=FALSE)

dev.off()



#---------------------------------------------------------------------
# Matrix Scatter Plots
#----------------------------------------------------------------------

## scatter plot for possible collinrarity
png("mt_splom_scatter_matrix.png",width=1000,height=800)
splom( data0[,c("X1","X2","muhat","Vhat")] , pscales = 0,main="Matrix Scatter Plot",
	col.main='red', font.main=4, xlab="")
dev.off()

par(mar = c(4,4,1,0))



png("mt_ggp_corr_plot.png",width=1000,height=800)

ggpairs(data0[,c("X1","X2","muhat","Vhat")],
	title = "Correlation Plot",
	upper = list ( 
		mapping = ggplot2::aes(size = 16), 
		color = 'red'
	),
	lower = list(
    	continuous = "smooth",
    	combo = "facetdensity",
    	mapping = ggplot2::aes(color = muhat)
  	),
	axisLabels='show')

dev.off()

#--------------------------------------------------------------------------
# X1,X2 - muhat,Vhat Plots
#--------------------------------------------------------------------------

png("mt_scatter_detailed_plot.png",width=1000,height=800)

par(mfrow = c(2,2))

## Or you can first create an initial plot of one line
##         and then iteratively add the lines
##
##   below is an example to plot X1 vs. muhat for different X2 values
##
flag <- which(data0$X2 == 0)
plot(data0$X1[flag], data0$muhat[flag], type="l", xlim=range(data0$X1), ylim=range(data0$muhat),
	xlab="X1", ylab="Mean(Y)",col="blue")
for (j in 1:40){
  flag <- which(data0$X2 == 0.1*j)
  lines(data0$X1[flag], data0$muhat[flag])
}
title(main="X1 Vs Mean(Y) - for different X2",
	col.main='red', font.main=4,outer=FALSE)


## Or you can first create an initial plot of one line
##         and then iteratively add the lines
##
##   below is an example to plot X2 vs. muhat for different X1 values
##
flag <- which(data0$X1 == 0)
plot(data0$X2[flag], data0$muhat[flag], type="l", xlim=range(data0$X2), ylim=range(data0$muhat),
	xlab="X2", ylab="Mean(Y)",col="blue")
for (j in 1:70){
  flag <- which(data0$X1 == 0.1*j)
  lines(data0$X2[flag], data0$muhat[flag])
}
title(main="X2 Vs Mean(Y) - for different X1",
	col.main='red', font.main=4,outer=FALSE)



# variance vs X1 for different values of X2
flag <- which(data0$X2 == 0)
plot(data0$X1[flag], data0$Vhat[flag], type="l", xlim=range(data0$X1), ylim=range(data0$Vhat),
	xlab="X1", ylab="Variance(Y)",col="red")
for (j in 1:40){
  flag <- which(data0$X2 == 0.1*j)
  lines(data0$X1[flag], data0$Vhat[flag])
}
title(main="X1 Vs Variance(Y)- for different X2",
	col.main='red', font.main=4,outer=FALSE)

# variance vs X2 for different values of X1
flag <- which(data0$X1 == 0)
plot(data0$X2[flag], data0$Vhat[flag], type="l", xlim=range(data0$X2), ylim=range(data0$Vhat),
	xlab="X2", ylab="Variance(Y)",col="red")
for (j in 1:70){
  flag <- which(data0$X1 == 0.1*j)
  lines(data0$X2[flag], data0$Vhat[flag])
}
title(main="X2 Vs Variance(Y) - for different X1",
	col.main='red', font.main=4,outer=FALSE)

dev.off()



#-----------------------------------------------------------
# X1 vs X2 scatter plots for mean and variance
#------------------------------------------------------------

png("mt_matplot_x1_x2_mu_v_qt.png",width=1000,height=800)

m <- matrix(c(1,2,3,4,5),nrow = 5,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(.05,0.375,.05,0.375,.05))

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main="X1 X2 Vs Mean(Y)/Variance(Y) Quartile Bands", col.main='red', font.main=4 ,outer=TRUE)

par(mar = c(4,4,1,0))
plot(data0$X1,data0$X2,type='p', pch=21, col=as.numeric(muqt), 
	bg=as.numeric(muqt),
	xlab="X1", ylab="X2")
title(main="Mean(Y) Quartiles", col.main='red', font.main=4 ,outer=FALSE)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

par(mar = c(4,4,1,0))
plot(data0$X1,data0$X2,type='p', pch=21, col=as.numeric(Vqt), 
	bg=as.numeric(Vqt),
	xlab="X1", ylab="X2")
title(main="Variance(Y) Quartiles", col.main='red', font.main=4,outer=FALSE)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="bottom", inset=0, paste("Quartile",levels(as.factor(Vqt))),
	 fill=levels(as.factor(muqt)) ,cex=1, horiz = TRUE)

dev.off()




#--------------------------------------------------------------------
# Muhat to variance plot, is the variance uniform ( required for ols method)
# If variance is different heterodescacity - need to look at logit which
# does not assume homodescasity
# Not sure if nomality is met here too
#------------------------------------------------------------------
png("mt_matplot_mu_v_qt.png",width=1000,height=800)

m <- matrix(c(1,2,3),nrow = 3,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))plot

layout(mat = m,heights = c(.05,0.9,.05))

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main="Mean(Y) Vs Variance(Y)", col.main='red', font.main=4 ,outer=TRUE)

par(mar = c(4,4,1,0))
plot(data0$muhat,data0$Vhat,type='p', pch=21, col=as.numeric(Vqt), 
	bg=as.numeric(Vqt),
	xlab="Mean(Y) - muhat", ylab="Variance(Y) - Vhat")
title(main="Mean(Y) Vs Variance(Y) ", col.main='red', font.main=4 ,outer=FALSE)


par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="bottom", inset=0, paste("Variance(Y) Quartile",levels(as.factor(Vqt))),
	 fill=levels(as.factor(muqt)) ,cex=1, horiz = TRUE)


dev.off()

	

#----------------------------------------------------------------
#  Fitting Models - Training and Cross Validation
#-----------------------------------------------------------------

# Keep track of MSE errors for different models
muhat.models <- c ('Logistic Regression','Logistic Regression(Poly)',
			'Linear Regression','Linear Regression(Poly)','LASSO')

vhat.models <- c ('Logistic Regression','Logistic Regression(Poly)',
			'Linear Regression','Linear Regression(Poly)','LASSO')

error.type <- c('TRAIN_MSE','TEST_MSE')

mse.models <- matrix(NA,length(muhat.models),length(error.type))
rownames(mse.models) <- muhat.models
colnames(mse.models) <- error.type

v.models <- matrix(NA,length(vhat.models),length(error.type))
rownames(v.models) <- vhat.models
colnames(v.models) <- error.type

# Split the training data into training and test
validation.test.percent = 5
test.count <- round(dim(data0)[1] * (validation.test.percent/100))
test.rows <- sort(sample(1:dim(data0)[1],test.count,replace=FALSE))

data0.test <- data0[test.rows,]
data0.train <- data0[-test.rows,]

muqt.test <- muqt[test.rows]
Vqt.test <- Vqt[test.rows]
muqt.train <- muqt[-test.rows]
Vqt.train <- Vqt[-test.rows]


# 10 Fold Cross Validation - Folds
#
nfolds <- 10
folds <- sample(1:nfolds,length(data0.train$X1),replace=TRUE)


#-------------------------------------------------------------------
#  Fitting muhat
#
#-------------------------------------------------------------------



#--------------------------------------------------------
# muhat - Standard Logistic regression
#-------------------------------------------------------
# train using all train data

# train using all train data for minimum poly
lr.mu.fit <- glm(muhat~X1*X2,
			data=data0.train,
			family=binomial(logit))

train.residuals <- data0.train$muhat-lr.mu.fit$fitted
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lr.mu.fit,data0.test[,c("X1","X2")],type='response')
test.residuals <- data0.test$muhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

mse.models["Logistic Regression",1] <- mean(train.residuals^2)
mse.models["Logistic Regression",2] <- mean(test.residuals^2)

summary(lr.mu.fit)

# chi sq test
chi.lr.muhat <- chisq.test(data0.train$muhat, lr.mu.fit$fitted)
print(chi.lr.muhat)

#-----------------------------------------------------------
# muhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_mean_trg_lr",
	data0.train,
	lr.mu.fit$residuals,
	train.stdres,
	muqt.train,
	data0.train$muhat,
	lr.mu.fit$fitted,
	"Mean(Y)",
	"Mean(Y) Quartile",
	"Logistic Regression"
)




#--------------------------------------------------------
# muhat - Poly Logistic regression
#-------------------------------------------------------
# train using all train data
# use 10 fold cross validation to choose the best poly term
poly.mse <- matrix(NA,poly.x1.max,poly.x2.max)
rownames(poly.mse) <- paste('X1_',c(1:6),sep="")
colnames(poly.mse) <- paste('X2_',c(1:6),sep="")

for ( p.x1 in 1:poly.x1.max ) {

	for ( p.x2 in 1:poly.x2.max ) {

		pred.mse <- matrix(NA,nfolds,1)

		for ( i in 1:nfolds) {
	
			lrp.mu.fit <- glm(muhat~poly(X1,p.x1)*poly(X2,p.x2),data=data0.train[folds!=i,],
							family=binomial(logit))

			lrp.pred <- predict(lrp.mu.fit,data0.train[folds==i,c("X1","X2")],type='response')

			pred.mse[i,1] <- mean((lrp.pred-data0.train$muhat[folds==i])^2)
		}

		poly.mse[p.x1,p.x2] <- apply(pred.mse,2,mean)
	}
}

# minimum poly with complexity factored in by muliplying with log of poly+1
poly.lrp.min <- arrayInd(which.min(poly.mse),dim(poly.mse))

poly.lrp.muhat.min.x1 <- poly.lrp.min[1]
poly.lrp.muhat.min.x2 <- poly.lrp.min[2]

# train using all train data for minimum poly
lrp.mu.fit <- glm(muhat~poly(X1,poly.lrp.muhat.min.x1)*poly(X2,poly.lrp.muhat.min.x2),
							data=data0.train,
							family=binomial(logit))

train.residuals <- data0.train$muhat-lrp.mu.fit$fitted
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lrp.mu.fit,data0.test[,c("X1","X2")],type='response')
test.residuals <- data0.test$muhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

mse.models["Logistic Regression(Poly)",1] <- mean(train.residuals^2)
mse.models["Logistic Regression(Poly)",2] <- mean(test.residuals^2)

summary(lrp.mu.fit)

# chi sq test
chi.lrp.muhat <- chisq.test(data0.train$muhat, lrp.mu.fit$fitted)
print(chi.lrp.muhat)

#-----------------------------------------------------------
# muhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_mean_trg_lr_poly",
	data0.train,
	lrp.mu.fit$residuals,
	train.stdres,
	muqt.train,
	data0.train$muhat,
	lrp.mu.fit$fitted,
	"Mean(Y)",
	"Mean(Y) Quartile",
	"Logistic Regression(Poly)"
)


#--------------------------------------------------------------------
# Muhat to variance plot, is the variance uniform ( required for ols method)
# If variance is different heterodescacity - need to look at logit which
# does not assume homodescasity
# Not sure if nomality is met here too
#------------------------------------------------------------------

rocol.mse <- arrayInd(c(1:length(c(poly.mse))),dim(poly.mse))
cv.labels <- paste(rownames(poly.mse)[rocol.mse[,1]],colnames(poly.mse)[rocol.mse[,2]])

png("mt_cvplot_mu_lg_poly.png",width=1000,height=800)

m <- matrix(c(1,2,3,3),nrow = 4,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(.05,0.80,.05,0.1))

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main="Logistic Regression Mean(Y)- Cross Validation MSE plot Vs Poly Terms", 
##	col.main='red', font.main=4 ,outer=TRUE)

par(mar = c(4,4,1,0))

plot(c(1:length(c(poly.mse))),c(poly.mse),xaxt='n',type='p',pch=21,
	xlab="Poly Terms", ylab="Cross Validation MSE - muhat", bg='blue')

lo <- loess(c(poly.mse)~c(1:length(c(poly.mse))))

lines(predict(lo), col='red', lwd=2)

axis(1,at=1:36,labels=cv.labels,cex=.5)


par(mar = c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="bottom", inset=0, 
	c("MSE from Logistic Regression Cross for Mean(Y)","Smoothing Line for MSE"),
	fill=c("blue","red") ,
	cex=1, horiz = FALSE)


dev.off()













#--------------------------------------------------------
# muhat - Standard linear regression
#-------------------------------------------------------
# train using all train data
lm.mu.fit <- lm(muhat~X1+X2,data=data0.train)

train.residuals <- data0.train$muhat-lm.mu.fit$fitted.values
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lm.mu.fit,data0.test[,c("X1","X2")])
test.residuals <- data0.test$muhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

mse.models['Linear Regression',1] <- mean(train.residuals^2)
mse.models['Linear Regression',2] <- mean(test.residuals^2)

summary(lm.mu.fit)

# chi sq test
chi.lm.muhat <- chisq.test(data0.train$muhat, lm.mu.fit$fitted.values)
print(chi.lm.muhat)

#-----------------------------------------------------------
# muhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_mean_trg_slg",
	data0.train,
	lm.mu.fit$residuals,
	train.stdres,
	muqt.train,
	data0.train$muhat,
	lm.mu.fit$fitted.values,
	"Mean(Y)",
	"Mean(Y) Quartile",
	"Linear Regression"
)


#--------------------------------------------------------
# muhat - linear regression  - Poly
#-------------------------------------------------------
# use 10 fold cross validation to choose the best poly term
poly.mse <- matrix(NA,poly.x1.max,poly.x2.max)

for ( p.x1 in 1:poly.x1.max ) {

	for ( p.x2 in 1:poly.x2.max ) {

		pred.mse <- matrix(NA,nfolds,1)

		for ( i in 1:nfolds) {
	
			lmp.mu.fit <- lm(muhat~poly(X1,p.x1)+poly(X2,p.x2),data=data0.train[folds!=i,])
		
			lmp.pred <- predict(lmp.mu.fit,data0.train[folds==i,c("X1","X2")])

			pred.mse[i,1] <- mean((lmp.pred-data0.train$muhat[folds==i])^2)
		}

		poly.mse[p.x1,p.x2] <- apply(pred.mse,2,mean)
	}
}

# minimum poly with complexity factored in by muliplying with log of poly+1
poly.min <- arrayInd(which.min(poly.mse),dim(poly.mse))

poly.muhat.min.x1 <- poly.min[1]
poly.muhat.min.x2 <- poly.min[2]

# train using all train data for minimum poly
lmp.mu.fit <- lm(muhat~poly(X1,poly.muhat.min.x1)+poly(X2,poly.muhat.min.x2),data=data0.train)

train.residuals <- data0.train$muhat-lmp.mu.fit$fitted.values
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lmp.mu.fit,data0.test[,c("X1","X2")])
test.residuals <- data0.test$muhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

mse.models['Linear Regression(Poly)',1] <- mean(train.residuals^2)
mse.models['Linear Regression(Poly)',2] <- mean(test.residuals^2)

summary(lmp.mu.fit)

# chi squared test for test data
chi.poly.muhat <- chisq.test(data0.train$muhat, lmp.mu.fit$fitted.values)
print(chi.poly.muhat)

#-----------------------------------------------------------
# muhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_mean_trg_lg_poly",
	data0.train,
	lmp.mu.fit$residuals,
	train.stdres,
	muqt.train,
	data0.train$muhat,
	lmp.mu.fit$fitted.values,
	"Mean(Y)",
	"Mean(Y) Quartile",
	"Linear Regression Polynomial"
)


#-----------------------------------------------------------
# muhat - Lasso
#-----------------------------------------------------------
# cross validation to choose lambda for lasso on the bets poly model
#
# 

lasso.lambda.grid <- 10^ seq (10,-6, length =1000)

poly.mse <- matrix(NA,poly.x1.max,poly.x2.max)

for ( p.x1 in 1:poly.x1.max ) {

	lasso.muhat.cols <- c(1:(2+p.x1-1))

	for ( p.x2 in 1:poly.x2.max ) {
	
		if ( p.x2 > 1 ) {
			lasso.muhat.cols <- c( lasso.muhat.cols,(2+poly.x1.max):(2+poly.x1.max+
							p.x2-2))
		}

		pred.mse <- matrix(NA,nfolds,1)

		for ( i in 1:nfolds) {
	
			lasso.mu.fit <- cv.glmnet(as.matrix(data0.train[folds!=i,lasso.muhat.cols]),
							data0.train$muhat[folds!=i],
							alpha=1, lambda=lasso.lambda.grid)

			lasso.muhat.bestlam <- lasso.mu.fit$lambda.min

			lasso.pred=predict(lasso.mu.fit, s=lasso.muhat.bestlam , 
						newx=as.matrix(data0.train[folds==i,lasso.muhat.cols]))

			pred.mse[i,1] <- mean((lasso.pred-data0.train$muhat[folds==i])^2)
		}

		poly.mse[p.x1,p.x2] <- apply(pred.mse,2,mean)
	}
}

# minimum poly with complexity factored in by muliplying with log of poly+1
poly.min <- arrayInd(which.min(poly.mse),dim(poly.mse))

poly.lasso.muhat.min.x1 <- poly.min[1]
poly.lasso.muhat.min.x2 <- poly.min[2]



# CV the best poly term with lasso on the whole data to get the best lambda for it, 
# using a wider range of lambda here
lasso.lambda.grid <- 10^ seq (10,-6, length =10000)
lasso.muhat.cols <- c(1:(2+poly.lasso.muhat.min.x1-1))
if ( poly.lasso.muhat.min.x2 > 1 ) {
	lasso.muhat.cols <- c( lasso.muhat.cols,(2+poly.x1.max):(2+poly.x1.max+poly.lasso.muhat.min.x2-2))
}

# cross validate to get the best lambda
lasso.mu.fit <- cv.glmnet(as.matrix(data0.train[,lasso.muhat.cols]),data0.train$muhat,
			alpha=1, lambda=lasso.lambda.grid)

lasso.muhat.bestlam <- lasso.mu.fit$lambda.min
coef(lasso.mu.fit, s = "lambda.min")

train.pred=predict(lasso.mu.fit, s=lasso.muhat.bestlam , newx=as.matrix(data0.train[,lasso.muhat.cols]))
train.residuals <- data0.train$muhat-train.pred
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred=predict(lasso.mu.fit, s=lasso.muhat.bestlam , newx=as.matrix(data0.test[,lasso.muhat.cols]))
test.residuals <- data0.test$muhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

mse.models['LASSO',1] <- mean(train.residuals^2)
mse.models['LASSO',2] <- mean(test.residuals^2)

summary(lasso.mu.fit)

tss <- sum((data0.train$muhat-mean(data0.train$muhat))^2)
rss <- sum(train.residuals^2)
lasso.muhat.R.squared <- (tss-rss)/tss
print(lasso.muhat.R.squared)


# chi squared test for test data
chi.lasso.muhat <- chisq.test(data0.train$muhat, train.pred)
print(chi.lasso.muhat)

#----------------------------------------------------------
# muhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
png("lasso_cv_plot_mean",width=1000,height=800)
plot(lasso.mu.fit)
title(main="Lasso CV Plot - Mean(Y)",
	col.main='red', font.main=4,outer=FALSE)
dev.off()


plotxy.residuals("mt_rse_plot_mean_trg_lasso_poly",
	data0.train,
	lasso.mu.fit$residuals,
	train.stdres,
	muqt.train,
	data0.train$muhat,
	train.pred,
	"Mean(Y)",
	"Mean(Y) Quartile",
	"Lasso Polynomial"
)



#-------------------------------------------------------------------
#  Fit models for Vhat
#
#-------------------------------------------------------------------



#--------------------------------------------------------
# vhat - Standard Logistic regression
#-------------------------------------------------------
# train using all train data

# train using all train data for minimum poly
lr.v.fit <- glm(Vhat~X1*X2,
			data=data0.train,
			family=binomial(logit))

train.residuals <- data0.train$Vhat-lr.v.fit$fitted
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lr.v.fit,data0.test[,c("X1","X2")],type='response')
test.residuals <- data0.test$muhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

v.models["Logistic Regression",1] <- mean(train.residuals^2)
v.models["Logistic Regression",2] <- mean(test.residuals^2)

summary(lr.v.fit)

# chi sq test
chi.lr.vhat <- chisq.test(data0.train$Vhat, lr.v.fit$fitted)
print(chi.lr.vhat)

#-----------------------------------------------------------
# vhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_v_trg_lr",
	data0.train,
	lr.v.fit$residuals,
	train.stdres,
	Vqt.train,
	data0.train$Vhat,
	lr.v.fit$fitted,
	"Variance(Y)",
	"Variance(Y) Quartile",
	"Logistic Regression"
)




#--------------------------------------------------------
# vhat - Poly Logistic regression
#-------------------------------------------------------
# train using all train data
# use 10 fold cross validation to choose the best poly term
poly.mse <- matrix(NA,poly.x1.max,poly.x2.max)
rownames(poly.mse) <- paste('X1_',c(1:6),sep="")
colnames(poly.mse) <- paste('X2_',c(1:6),sep="")

for ( p.x1 in 1:poly.x1.max ) {

	for ( p.x2 in 1:poly.x2.max ) {

		pred.mse <- matrix(NA,nfolds,1)

		for ( i in 1:nfolds) {
	
			lrp.v.fit <- glm(Vhat~poly(X1,p.x1)*poly(X2,p.x2),data=data0.train[folds!=i,],
							family=binomial(logit))

			lrp.pred <- predict(lrp.v.fit,data0.train[folds==i,c("X1","X2")],type='response')

			pred.mse[i,1] <- mean((lrp.pred-data0.train$Vhat[folds==i])^2)
		}

		poly.mse[p.x1,p.x2] <- apply(pred.mse,2,mean)
	}
}

# minimum poly with complexity factored in by muliplying with log of poly+1
poly.lrp.min <- arrayInd(which.min(poly.mse),dim(poly.mse))

poly.lrp.vhat.min.x1 <- poly.lrp.min[1]
poly.lrp.vhat.min.x2 <- poly.lrp.min[2]

# train using all train data for minimum poly
lrp.v.fit <- glm(Vhat~poly(X1,poly.lrp.vhat.min.x1)*poly(X2,poly.lrp.vhat.min.x2),
							data=data0.train,
							family=binomial(logit))

train.residuals <- data0.train$Vhat-lrp.v.fit$fitted
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lrp.v.fit,data0.test[,c("X1","X2")],type='response')
test.residuals <- data0.test$Vhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

v.models["Logistic Regression(Poly)",1] <- mean(train.residuals^2)
v.models["Logistic Regression(Poly)",2] <- mean(test.residuals^2)

summary(lrp.v.fit)

# chi sq test
chi.lrp.vhat <- chisq.test(data0.train$Vhat, lrp.v.fit$fitted)
print(chi.lrp.vhat)

#-----------------------------------------------------------
# vhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_v_trg_lr_poly",
	data0.train,
	lrp.v.fit$residuals,
	train.stdres,
	Vqt.train,
	data0.train$Vhat,
	lrp.v.fit$fitted,
	"Variance(Y)",
	"Variance(Y) Quartile",
	"Logistic Regression(Poly)"
)


#--------------------------------------------------------------------
# Muhat to variance plot, is the variance uniform ( required for ols method)
# If variance is different heterodescacity - need to look at logit which
# does not assume homodescasity
# Not sure if nomality is met here too
#------------------------------------------------------------------

rocol.mse <- arrayInd(c(1:length(c(poly.mse))),dim(poly.mse))
cv.labels <- paste(rownames(poly.mse)[rocol.mse[,1]],colnames(poly.mse)[rocol.mse[,2]])

png("mt_cvplot_v_lg_poly.png",width=1000,height=800)

m <- matrix(c(1,2,3,3),nrow = 4,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(.05,0.80,.05,0.1))

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

par(mar = c(4,4,1,0))

plot(c(1:length(c(poly.mse))),c(poly.mse),xaxt='n',type='p',pch=21,
	xlab="Poly Terms", ylab="Cross Validation MSE - Vhat- Variance(Y)", bg='blue')

lo <- loess(c(poly.mse)~c(1:length(c(poly.mse))))

lines(predict(lo), col='red', lwd=2)

axis(1,at=1:36,labels=cv.labels,cex=.5)


par(mar = c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="bottom", inset=0, 
	c("MSE from Logistic Regression Cross for Variance(Y)","Smoothing Line for MSE"),
	fill=c("blue","red") ,
	cex=1, horiz = FALSE)


dev.off()





































#--------------------------------------------------------
# Vhat  - Standard linear regression
#-------------------------------------------------------
# train using all train data
lm.v.fit <- lm(Vhat~X1+X2,data=data0.train)

train.residuals <- data0.train$Vhat-lm.v.fit$fitted.values
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lm.v.fit,data0.test[,c("X1","X2")])
test.residuals <- data0.test$Vhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

v.models['Linear Regression',1] <- mean(train.residuals^2)
v.models['Linear Regression',2] <- mean(test.residuals^2)

summary(lm.v.fit)


# chi squared test for test data
chi.lm.vhat <- chisq.test(data0.train$Vhat, lm.v.fit$fitted.values)
print(chi.lm.vhat)

#-----------------------------------------------------------
# Plot residuals vs X1,X2 by Variance
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_var_v_trg_slg",
	data0.train,
	lm.v.fit$residuals,
	train.stdres,
	Vqt.train,
	data0.train$Vhat,
	lm.v.fit$fitted.values,
	"Variance(Y)",
	"Variance(Y) Quartile",
	"Linear Regression"
)



#--------------------------------------------------------
# Vhat - linear regression  - Poly
#-------------------------------------------------------
# use 10 fold cross validation to choose the best poly term
poly.mse <- matrix(NA,poly.x1.max,poly.x2.max)

for ( p.x1 in 1:poly.x1.max ) {

	for ( p.x2 in 1:poly.x2.max ) {

		pred.mse <- matrix(NA,nfolds,1)

		for ( i in 1:nfolds) {
	
			lmp.v.fit <- lm(Vhat~poly(X1,p.x1)+poly(X2,p.x2),data=data0.train[folds!=i,])
		
			lmp.pred <- predict(lmp.v.fit,data0.train[folds==i,c("X1","X2")])

			pred.mse[i,1] <- mean((lmp.pred-data0.train$Vhat[folds==i])^2)
		}

		poly.mse[p.x1,p.x2] <- apply(pred.mse,2,mean)
	}
}

# minimum poly with complexity factored in by muliplying with log of poly+1
poly.min <- arrayInd(which.min(poly.mse),dim(poly.mse))

poly.Vhat.min.x1 <- poly.min[1]
poly.Vhat.min.x2 <- poly.min[2]

# train using all train data for minimum poly
lmp.v.fit <- lm(Vhat~poly(X1,poly.Vhat.min.x1)+poly(X2,poly.Vhat.min.x2),data=data0.train)

train.residuals <- data0.train$Vhat-lmp.v.fit$fitted.values
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred <- predict(lmp.v.fit,data0.test[,c("X1","X2")])
test.residuals <- data0.test$Vhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

v.models['Linear Regression(Poly)',1] <- mean(train.residuals^2)
v.models['Linear Regression(Poly)',2] <- mean(test.residuals^2)

summary(lmp.v.fit)

# chi squared test for test data
chi.poly.vhat <- chisq.test(data0.train$Vhat, lmp.v.fit$fitted.values)
print(chi.poly.vhat)


#-----------------------------------------------------------
# Vhat - Plot residuals vs X1,X2 by Variance
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_var_v_trg_lg_poly",
	data0.train,
	lmp.v.fit$residuals,
	train.stdres,
	Vqt.train,
	data0.train$Vhat,
	lmp.v.fit$fitted.values,
	"Variance(Y)",
	"Variance(Y) Quartile",
	"Linear Regression Polynomial - Variance(Y)"
)



#-----------------------------------------------------------
# Vhat - Lasso
#-----------------------------------------------------------
# cross validation to choose lambda for lasso on the bets poly model
#
# 

lasso.lambda.grid <- 10^ seq (10,-6, length =1000)

poly.mse <- matrix(NA,poly.x1.max,poly.x2.max)

for ( p.x1 in 1:poly.x1.max ) {

	lasso.Vhat.cols <- c(1:(2+p.x1-1))

	for ( p.x2 in 1:poly.x2.max ) {
	
		if ( p.x2 > 1 ) {
			lasso.Vhat.cols <- c( lasso.Vhat.cols,(2+poly.x1.max):(2+poly.x1.max+
							p.x2-2))
		}

		pred.mse <- matrix(NA,nfolds,1)

		for ( i in 1:nfolds) {
	
			lasso.v.fit <- cv.glmnet(as.matrix(data0.train[folds!=i,lasso.Vhat.cols]),
							data0.train$Vhat[folds!=i],
							alpha=1, lambda=lasso.lambda.grid)

			lasso.Vhat.bestlam <- lasso.v.fit$lambda.min

			lasso.pred=predict(lasso.v.fit, s=lasso.Vhat.bestlam , 
						newx=as.matrix(data0.train[folds==i,lasso.Vhat.cols]))

			pred.mse[i,1] <- mean((lasso.pred-data0.train$Vhat[folds==i])^2)
		}

		poly.mse[p.x1,p.x2] <- apply(pred.mse,2,mean)
	}
}

# minimum poly with complexity factored in by muliplying with log of poly+1
poly.min <- arrayInd(which.min(poly.mse),dim(poly.mse))

poly.lasso.Vhat.min.x1 <- poly.min[1]
poly.lasso.Vhat.min.x2 <- poly.min[2]



# CV the best poly term with lasso on the whole data to get the best lambda for it, 
# using a wider range of lambda here
lasso.lambda.grid <- 10^ seq (10,-6, length =10000)
lasso.Vhat.cols <- c(1:(2+poly.lasso.Vhat.min.x1-1))
if ( poly.lasso.Vhat.min.x2 > 1 ) {
	lasso.Vhat.cols <- c( lasso.Vhat.cols,(2+poly.x1.max):(2+poly.x1.max+poly.lasso.Vhat.min.x2-2))
}

# cross validate to get the best lambda
lasso.v.fit <- cv.glmnet(as.matrix(data0.train[,lasso.Vhat.cols]),data0.train$Vhat,
			alpha=1, lambda=lasso.lambda.grid)

lasso.Vhat.bestlam <- lasso.v.fit$lambda.min
coef(lasso.v.fit, s = "lambda.min")

train.pred=predict(lasso.v.fit, s=lasso.Vhat.bestlam , newx=as.matrix(data0.train[,lasso.Vhat.cols]))
train.residuals <- data0.train$Vhat-train.pred
train.sigma <- sqrt(sum(train.residuals^2)/(dim(data0.train)[1]-2-1))
train.stdres <- train.residuals / train.sigma

test.pred=predict(lasso.v.fit, s=lasso.Vhat.bestlam , newx=as.matrix(data0.test[,lasso.Vhat.cols]))
test.residuals <- data0.test$Vhat - test.pred
test.sigma <- sqrt(sum(test.residuals^2)/(dim(data0.test)[1]-2-1))
test.stdres <- test.residuals / test.sigma

v.models['LASSO',1] <- mean(train.residuals^2)
v.models['LASSO',2] <- mean(test.residuals^2)

summary(lasso.v.fit)

tss <- sum((data0.train$Vhat-mean(data0.train$Vhat))^2)
rss <- sum(train.residuals^2)
lasso.Vhat.R.squared <- (tss-rss)/tss
print(lasso.Vhat.R.squared)


# chi squared test for test data
chi.lasso.vhat <- chisq.test(data0.train$Vhat, train.pred)
print(chi.lasso.vhat)

#----------------------------------------------------------
# Vhat - Plot residuals vs X1,X2 by mean
#----------------------------------------------------------
png("lasso_cv_plot_v.png",width=1000,height=800)
plot(lasso.v.fit)
title(main="Lasso CV Plot - Variance(Y)",
	col.main='red', font.main=3)
dev.off()


#-----------------------------------------------------------
# Vhat Plot residuals vs X1,X2 by Variance
#----------------------------------------------------------
plotxy.residuals("mt_rse_plot_var_v_trg_lasso_poly",
	data0.train,
	lasso.v.fit$residuals,
	train.stdres,
	Vqt.train,
	data0.train$Vhat,
	train.pred,
	"Variance(Y)",
	"Variance(Y) Quartile",
	"Lasso - Variance(Y)"
)




#-----------------------------------------------------------------------
# Plot the MSE test error for different models
#-----------------------------------------------------------------------
print(mse.models)

png("mt_train_test_mse_muhat_model_plot.png",width=1000,height=800)

matplot(mse.models, type="b", lty=1,xlab="Model", ylab="MSE(Mean(Y))",
	col=seq_len(ncol(mse.models)),
	cex.axis=1,xaxt="n")
axis(1, at = 1:length(muhat.models), labels = paste(muhat.models), cex.axis = .5)
title(main="Model Vs MSE(Mean(Y)) Training/Test",col.main='red', font.main=4)
legend("topright", colnames(mse.models),col=seq_len(ncol(mse.models)),cex=0.8,
	fill=seq_len(ncol(mse.models)))

dev.off()



print(v.models)

png("mt_train_test_mse_v_model_plot.png",width=1000,height=800)

matplot(v.models, type="b", lty=1,xlab="Model", ylab="MSE(Variance(Y))",
	col=seq_len(ncol(v.models)),
	cex.axis=1,xaxt="n")
axis(1, at = 1:length(vhat.models), labels = paste(vhat.models), cex.axis = .5)
title(main="Model Vs MSE(Variance(Y)) Training/Test",col.main='red', font.main=4)
legend("topright", colnames(v.models),col=seq_len(ncol(v.models)),cex=0.8,
	fill=seq_len(ncol(v.models)))

dev.off()




#----------------------------------------------------------------------
# Do a Bootstrap test to pick the best model statistically
#----------------------------------------------------------------------
### number of loops
B <- 100

### Final TE values for Mean(Y)
TEALL <- NULL

### Final TE values for Variance(y)
VEALL <- NULL

# create a matrix to hold the test MSE of the different models for each of the B cycles
TEALL <- matrix(NA,B,length(muhat.models),dimnames=list(1:B,muhat.models))
VEALL <- matrix(NA,B,length(vhat.models),dimnames=list(1:B,vhat.models));

# Split the training data into training and test
validation.test.percent <- 50

for (b in 1:B){

	### randomly select 10% observations as testing data in each loop
	test.count <- round(dim(data0)[1] * (validation.test.percent/100))
	test.rows <- sort(sample(1:dim(data0)[1],test.count,replace=TRUE))

	data0.test <- data0[test.rows,]
	data0.train <- data0[-test.rows,]

	muqt.test <- muqt[test.rows]
	Vqt.test <- Vqt[test.rows]
	muqt.train <- muqt[-test.rows]
	Vqt.train <- Vqt[-test.rows]



	#---------------------------------------
	#  Logistic Regression
	#---------------------------------------
	# muhat train using all train data
	lr.mu.fit <- glm(muhat~X1*X2,
				data=data0.train,
				family=binomial(logit))

	test.pred <- predict(lr.mu.fit,data0.test[,c("X1","X2")],type='response')
	test.residuals <- data0.test$muhat - test.pred
	TEALL[b,"Logistic Regression"] <- mean(test.residuals^2)

	# vhat train using all train data
	lr.v.fit <- glm(Vhat~X1*X2,
				data=data0.train,
				family=binomial(logit))

	test.pred <- predict(lr.v.fit,data0.test[,c("X1","X2")],type='response')
	test.residuals <- data0.test$Vhat - test.pred
	VEALL[b,"Logistic Regression"] <- mean(test.residuals^2)


	#---------------------------------------
	#  Logistic Regression Poly
	#---------------------------------------
	# muhat train using all train data
	lrp.mu.fit <- glm(muhat~poly(X1,poly.lrp.muhat.min.x1)*poly(X2,poly.lrp.muhat.min.x2),
							data=data0.train,
							family=binomial(logit))

	test.pred <- predict(lrp.mu.fit,data0.test[,c("X1","X2")],type='response')
	test.residuals <- data0.test$muhat - test.pred
	TEALL[b,"Logistic Regression(Poly)"] <- mean(test.residuals^2)


	# vhat train using all train data
	lrp.v.fit <- glm(Vhat~poly(X1,poly.lrp.vhat.min.x1)*poly(X2,poly.lrp.vhat.min.x2),
							data=data0.train,
							family=binomial(logit))

	test.pred <- predict(lrp.v.fit,data0.test[,c("X1","X2")],type='response')
	test.residuals <- data0.test$Vhat - test.pred
	VEALL[b,"Logistic Regression(Poly)"] <- mean(test.residuals^2)


	#---------------------------------------
	#  Linear Regression 
	#---------------------------------------
	# muhat train using all train data
	lm.mu.fit <- lm(muhat~X1+X2,data=data0.train)
	test.pred <- predict(lm.mu.fit,data0.test[,c("X1","X2")])
	test.residuals <- data0.test$muhat - test.pred
	TEALL[b,"Linear Regression"] <- mean(test.residuals^2)

	# Vhat train using all train data
	lm.v.fit <- lm(Vhat~X1+X2,data=data0.train)
	test.pred <- predict(lm.v.fit,data0.test[,c("X1","X2")])
	test.residuals <- data0.test$Vhat - test.pred
	VEALL[b,"Linear Regression"] <- mean(test.residuals^2)



	#---------------------------------------
	#  Linear Regression Poly
	#---------------------------------------
	# muhat train using all train data
	lmp.mu.fit <- lm(muhat~poly(X1,poly.muhat.min.x1)+poly(X2,poly.muhat.min.x2),data=data0.train)
	test.pred <- predict(lmp.mu.fit,data0.test[,c("X1","X2")])
	test.residuals <- data0.test$muhat - test.pred
	TEALL[b,"Linear Regression(Poly)"] <- mean(test.residuals^2)

	# Vhat train using all train data
	lmp.v.fit <- lm(Vhat~poly(X1,poly.Vhat.min.x1)+poly(X2,poly.Vhat.min.x2),data=data0.train)
	test.pred <- predict(lmp.v.fit,data0.test[,c("X1","X2")])
	test.residuals <- data0.test$Vhat - test.pred
	VEALL[b,"Linear Regression(Poly)"] <- mean(test.residuals^2)


	#---------------------------------------
	#  Lasso
	#---------------------------------------
	# muhat
	lasso.mu.fit <- glmnet(as.matrix(data0.train[,lasso.muhat.cols]),data0.train$muhat,alpha=1, 
					lambda=lasso.muhat.bestlam)
	test.pred <- predict(lasso.mu.fit, s=lasso.muhat.bestlam , 
					newx=as.matrix(data0.test[,lasso.muhat.cols]))
	test.residuals <- data0.test$muhat - test.pred
	TEALL[b,"LASSO"] <- mean(test.residuals^2)

	# Vhat
	lasso.v.fit <- glmnet(as.matrix(data0.train[,lasso.Vhat.cols]),data0.train$Vhat,alpha=1, 
					lambda=lasso.Vhat.bestlam)
	test.pred <- predict(lasso.v.fit, s=lasso.Vhat.bestlam , 
				newx=as.matrix(data0.test[,lasso.Vhat.cols]))
	test.residuals <- data0.test$Vhat - test.pred
	VEALL[b,"LASSO"] <- mean(test.residuals^2)


}


#---------------------------------------------------------------------
#box plots of the Bootstrap B=100 run results
#----------------------------------------------------------------------

TM <- as.matrix(apply(TEALL, 2, mean))
TS <- as.matrix(apply(TEALL, 2, var))

VM <- as.matrix(apply(VEALL, 2, mean))
VS <- as.matrix(apply(VEALL, 2, var))

colnames(TM) <- "Average MSE(Mean(Y))"
colnames(TS) <- "Variance MSE(Mean(Y))"

colnames(VM) <- "Average MSE(Variance(Y))"
colnames(VS) <- "Variance MSE(Variance(Y))"

print(TM)
print(TS)

print(VM)
print(VS)


#---------------------------------------------------------------------
# PLOT
#---------------------------------------------------------------------
# muhat avergae MSE

png("mt_bootstrap_muhat_mse_var_plot.png",width=1000,height=800)

m <- matrix(c(1,2,3),nrow = 3,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(0.05,0.475,0.475))

par(mar = c(0,0,0,5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main='BootStrap B=100 Results - Mean(Y)', col.main='red', font.main=4, cex=1.5, outer=TRUE)

par(mar = c(4,4,1,0))

matplot(TM, type="b", lty=1, xlab="Model", ylab="Average MSE(Mean(Y))",
	col=seq_len(nrow(TM)),pch=1,xaxt="n")
axis(1, at = 1:length(muhat.models), labels = paste(muhat.models), cex.axis = 1)
title(main="Bootstrap Average MSE for Mean(Y)",
	col.main='red', font.main=4)
legend("topright", colnames(TM),col=seq_len(ncol(TM)),cex=0.8,fill=seq_len(ncol(TM)))


## Variance 
par(mar = c(4,4,1,0))

matplot(TS, type="b", lty=1, xlab="Model", ylab="Variance MSE(Mean(Y))",
	col=seq_len(nrow(TS)),pch=1,xaxt="n")
axis(1, at = 1:length(muhat.models), labels = paste(muhat.models), cex.axis = 1)
title(main="Bootstrap Variance MSE for Mean(Y)",
	col.main='red', font.main=4)
legend("topright", colnames(TS),col=seq_len(ncol(TS)),cex=0.8,fill=seq_len(ncol(TS)))

dev.off()



# Vhat avergae MSE

png("mt_bootstrap_vhat_mse_var_plot.png",width=1000,height=800)

m <- matrix(c(1,2,3),nrow = 3,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(0.05,0.475,0.475))

par(mar = c(0,0,0,5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main='BootStrap B=100 Results - Variance(Y)', col.main='red', font.main=4, cex=1.5, outer=TRUE)

par(mar = c(4,4,1,0))

matplot(VM, type="b", lty=1, xlab="Model", ylab="Average MSE(Variance(Y))",
	col=seq_len(nrow(VM)),pch=1,xaxt="n")
axis(1, at = 1:length(vhat.models), labels = paste(vhat.models), cex.axis = 1)
title(main="Bootstrap Average MSE for Variance(Y)",
	col.main='red', font.main=4)
legend("topright", colnames(VM),col=seq_len(ncol(VM)),cex=0.8,fill=seq_len(ncol(VM)))


## Variance 
par(mar = c(4,4,1,0))

matplot(VS, type="b", lty=1, xlab="Model", ylab="Variance MSE(Variance(Y))",
	col=seq_len(nrow(VS)),pch=1,xaxt="n")
axis(1, at = 1:length(vhat.models), labels = paste(vhat.models), cex.axis = 1)
title(main="Bootstrap Variance MSE for Variance(Y)",
	col.main='red', font.main=4)
legend("topright", colnames(VS),col=seq_len(ncol(VS)),cex=0.8,fill=seq_len(ncol(VS)))

dev.off()






##boxplot
png("mt_boxplot_bootstrap_muhat_results.png",width=1000,height=800)

m <- matrix(c(1,2,3,4),nrow = 4,ncol = 1,byrow = TRUE)

par(oma=c(1,1,1,1))

layout(mat = m,heights = c(0.05,0.45,0.05,0.45))

par(mar = c(0,0,0,5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##title(main='BootStrap B=100 Results', col.main='red', font.main=4, cex=1.5, outer=TRUE)

par(mar = c(4,4,1,0))
boxplot(TEALL,col=(c("gold","darkgreen")),main="Bootstrap MSE for Mean(Y) by Model", 
		xlab="Model",ylab="MSE Mean(Y)",
		col.main='red', font.main=4)

par(mar = c(0,0,0,5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

par(mar = c(4,4,1,0))
boxplot(VEALL,col=(c("gold","darkgreen")),main="Bootstrap MSE for Variance(Y) by Model", 
		xlab="Model",ylab="MSE Variance(Y)",
		col.main='red', font.main=4)

dev.off()



#-------------------------------------------------------------------------------
# T test and W test to pick the best method
#-------------------------------------------------------------------------------
# create a matrix to hold the test MSE of the different models for each of the B cycles

## mhat
best.mhat.method <- rownames(TM)[which.min(TM)]
print(best.mhat.method)

t.test.table.mhat <- matrix(NA,2,length(muhat.models)-1,dimnames=list(c('T-Test','Wilcox Test'), 
			muhat.models[muhat.models!=best.mhat.method]));

for ( m in muhat.models[muhat.models!=best.mhat.method] ) {

	t.test.table.mhat[1,m] <- t.test(TEALL[,best.mhat.method],TEALL[,m],paired=TRUE,
					conf.level=.95)$p.value
	t.test.table.mhat[2,m] <- wilcox.test(TEALL[,best.mhat.method],TEALL[,m],paired=TRUE,
				conf.level=.95)$p.value	
}

print(t.test.table.mhat)



## Vhat
best.Vhat.method <- rownames(VM)[which.min(VM)]
print(best.Vhat.method)

t.test.table.vhat <- matrix(NA,2,length(vhat.models)-1,dimnames=list(c('T-Test','Wilcox Test'), 
			vhat.models[vhat.models!=best.Vhat.method]));

for ( m in vhat.models[vhat.models!=best.Vhat.method] ) {

	t.test.table.vhat[1,m] <- t.test(VEALL[,best.Vhat.method],VEALL[,m],paired=TRUE,
					conf.level=.95)$p.value
	t.test.table.vhat[2,m] <- wilcox.test(VEALL[,best.Vhat.method],VEALL[,m],paired=TRUE,
				conf.level=.95)$p.value	
}

print(t.test.table.vhat)



#------------------------------------------------------------------------
# Pick the best model and train it on the whole training data
#------------------------------------------------------------------------'

X1_poly <- poly(distpredtest[,1],poly.x1.max,raw=TRUE)[,-1]
X2_poly <- poly(distpredtest[,2],poly.x2.max,raw=TRUE)[,-1]

distpredtest.datapoly <- data.frame(X1=distpredtest[,1],
	X2=distpredtest[,2],
	X1_poly,
	X2_poly)


# muhat train using all train data poly LM
#lmp.mu.fit <- lm(muhat~poly(X1,poly.muhat.min.x1)+poly(X2,poly.muhat.min.x2),data=data0)
#distpredtest.poly.mean <- predict(lmp.mu.fit,distpredtest.datapoly[,c("X1","X2")])


# Vhat train using all train data poly LM
#lmp.v.fit <- lm(Vhat~poly(X1,poly.Vhat.min.x1)+poly(X2,poly.Vhat.min.x2),data=data0)
#distpredtest.lmp.variance <- predict(lmp.v.fit,distpredtest.datapoly[,c("X1","X2")])



# T test and W test says Lasso is good too based on p-value
#lasso.v.fit <- glmnet(as.matrix(data0[,lasso.Vhat.cols]),data0$Vhat,alpha=1, 
#					lambda=lasso.Vhat.bestlam)

#distpredtest.lasso.variance <- predict(lasso.v.fit, s=lasso.Vhat.bestlam , 
#					newx=as.matrix(distpredtest.datapoly[,lasso.Vhat.cols]))


# Logistic regression muhat
# muhat train using all train data
	lrp.mu.fit <- glm(muhat~poly(X1,poly.lrp.muhat.min.x1)*poly(X2,poly.lrp.muhat.min.x1),
							data=data0,
							family=binomial(logit))

	distpredtest.mean <- predict(lrp.mu.fit,distpredtest.datapoly[,c("X1","X2")],
		type='response')

# vhat train using all train data
	lrp.v.fit <- glm(Vhat~poly(X1,poly.lrp.vhat.min.x1)*poly(X2,poly.lrp.vhat.min.x1),
							data=data0,
							family=binomial(logit))

	distpredtest.variance <- predict(lrp.v.fit,distpredtest.datapoly[,c("X1","X2")],
		type='response')





#-------------------------------------------------------------------------
# Write the results csv file
#-------------------------------------------------------------------------
resultdata <- data.frame(distpredtest,
	mean=format(round(distpredtest.mean,6),nsmall=6),
	variance=format(round(distpredtest.variance,6),nsmall=6))

write.csv(resultdata, file = "distpred_results.csv",quote=FALSE)

