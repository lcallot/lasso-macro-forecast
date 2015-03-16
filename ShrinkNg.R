library('grplasso')



# Increment the date based on monthly obs by 'inc' month
incdate	<-function(date,inc)
{
	m	<-date[2]
	y	<-date[1]

	mi<-m+inc
	yi<-y
	if(mi>12 )
	{
		yi	<-yi+mi%/%12
		mi	<-mi%%12
	}
	dateinc	<-c(yi,mi)

	return(dateinc)
}

# Lag selection for the OLS-VAR by BIC
VARBIC	<-function(sldata,maxlag,direct=FALSE,...)
{
	arglist<-list(...)
	if(direct){fchor<-arglist[[1]]}
	if(!direct)fchor<-1
	lBIC	<-rep(0,maxlag)

	for(l in 1:maxlag)
	{
		y	<-sldata
		nbrser	<-ncol(y)
		for(m in 1:l){y<-ts.intersect(y,lag(sldata,1-fchor-m))}
		mcres	<-mclapply(1:nbrser,lmRSS,y,cbind(1,y[,-c(1:nbrser)]))
		res	<-Reduce('+',mcres)
		lBIC[l]	<-log(res) + (nbrser)*(1+l*nbrser)*log(nrow(y))/nrow(y)	
	}
	return(which(lBIC==min(lBIC),arr.ind=TRUE))
}
	
#Used by VARBIC to get the RSS of the VAR
lmRSS<-function(i,y,x){return(sum((lm.fit(y=y[,i],x=x)$resid)^2))}



#	Lasso estimation using glmnet
mclasso	<-function(i,y,nbrser)
{
	lnet	<-glmnet(x=y[,-c(1:nbrser)],y=y[,i],family='gaussian',alpha=1,nlambda=50,standardize=TRUE,type.gaussian='covariance')
	return(lnet)
}

#	aLasso with glmnet
mcada	<-function(i,y,nbrser,wols)
{
	gnet	<-glmnet(x=y[,-c(1:nbrser)],y=y[,i],family='gaussian',alpha=1,nlambda=50,standardize=FALSE,type.gaussian='covariance',penalty.factor=wols[i,])
	return(gnet)
}

#	agLasso estimation
mcgrp	<-function(i,grpindex,nbrser,y,bols,gamma)
	{
	ols	<-	bols[,-1]
	gW	=	rep(1,length(grpindex))
	for(g in 1:max(grpindex)){
		w	<-	sqrt(sum(bols[i,which(grpindex==g)]^2))^(-gamma)
		gW[which(grpindex==g)]	<-w
	}
	X	<- t(t(y[,-c(1:nbrser)])/gW)

	lambda	<-	lambdamax(x=cbind(1,X),y=y[,i],model=LinReg(),index=c(NA,grpindex))
	gpada	<-	grplasso(x=cbind(1,X),y=y[,i],model=LinReg(),standardized=TRUE,index=c(NA,grpindex),lambda=lambda)
	return(gpada)
	}


#Group Lasso estimation
estgrp	<-function(grpindex,nbrser,y,gamma)
{
	bols		<-NULL
	bgp.all		<-NULL

	for(i in 1:nbrser){
		obs	<-	length(y)/nbrser
		olsb	<-	lm.fit(y=y[,i],x=cbind(1,y[,-c(1:nbrser)]))$coeff
		bols	<-	rbind(bols,olsb)
	}
	if(sum(is.na(bols))!=0)print(paste(sum(is.na(bols)),'Na found in init ols for grp, replaced by 0s',sep=''))
	bols[which(is.na(bols))]=0
	#Estimation	
	mcgrpest	<-	mclapply(1:nbrser,mcgrp,grpindex,nbrser,y,bols,gamma,mc.silent=TRUE)

	for(i in 1:nbrser)
	{
	gpada	<-mcgrpest[[i]]
	bgp	<-gpada$coef
	bgp.all	<-rbind(bgp.all,t(bgp))
	}
grp<-list('grp'=bgp.all,'y'=y)
return(grp)
}







# Estimate the model on the training sample
mcadaesti	<-function(whichvar,y,gamma,nbrser)
{
	bols	<-NULL
	ada.b	<-NULL
	las.b	<-NULL

	# Initial OLS and weights
	for(i in 1:nbrser)
	{	
		obs	<-	length(sldata)/nbrser
		olsb	<-	lm.fit(y=y[,i],x=cbind(1,y[,-c(1:nbrser)]))$coeff
		bols	<-	rbind(bols,olsb)
	}
	if(sum(is.na(bols))!=0)print(paste(sum(is.na(bols)),'Na found in init ols for ADA, replaced by 0s',sep=''))
	bols[which(is.na(bols))]=0
	lasso	<-mclapply(c(1:nbrser),mclasso,y,nbrser)
	gada	<-list()
	aplas <- matrix(0,nbrser,ncol(y)-nbrser+1)

	for(g in gamma){wols<-abs(bols[,-1])^(-g);glnet<-mclapply(c(1:nbrser),mcada,y,nbrser,wols);gada<-c(gada,glnet)}

	for(i in 1:nbrser)
	{
		ada.cst	<-	NULL
		ada.bet	<-	NULL
		ada.BIC	<-	NULL

		las.cst	<-	NULL
		las.bet	<-	NULL
		las.BIC	<-	NULL

	# Plain old lasso
		las.mod	<-lasso[[i]]
		las.pred<-predict(las.mod,y[,-c(1:nbrser)])
		las.res	<-las.pred-matrix(y[,i],nrow=nrow(y),ncol=ncol(las.pred))
	#Selection by BIC (lasso)
		las.RSS	<-colSums(las.res^2)
		las.BIC	<-c(las.BIC,log(las.RSS[]/las.mod$nobs) + las.mod$df[]*log(las.mod$nobs)/las.mod$nobs)
		las.cst	<-c(las.cst,las.mod$a0)
		las.bet	<-cbind(las.bet,as.matrix(las.mod$beta))

		for(g in gamma){
	# Ada lasso		
			ada.mod	<-gada[[i]]
			ada.pred<-predict(ada.mod,y[,-c(1:nbrser)])
			ada.res	<-ada.pred-matrix(y[,i],nrow=nrow(y),ncol=ncol(ada.pred))
	#Selection by BIC (ada)
			ada.RSS	<-colSums(ada.res^2)
			ada.BIC	<-c(ada.BIC,log(ada.RSS[]/ada.mod$nobs) + ada.mod$df[]*log(ada.mod$nobs)/ada.mod$nobs)
			ada.cst	<-c(ada.cst,ada.mod$a0)
			ada.bet	<-cbind(ada.bet,as.matrix(ada.mod$beta))
			}
		ada.bmin	<-ada.bet[,which.min(ada.BIC)]	
		ada.b		<-rbind(ada.b,c(ada.cst[which.min(ada.BIC)],ada.bet[,which.min(ada.BIC)]))

		las.bmin	<-las.bet[,which.min(las.BIC)]	
		las.b		<-rbind(las.b,c(las.cst[which.min(las.BIC)],las.bet[,which.min(las.BIC)]))
# estimation of the ada post lasso
		if(max(abs(las.bmin))>0) #checking that not all variables were excluded
		{
			sel	<-c(rep(0,nbrser),las.bmin[-1])
			x	<-matrix(y[,sel!=0],nrow=nrow(y))
			wlas	<- abs((las.bmin[-1])[las.bmin[-1]!=0])^(-1)
			apl	<-glmnet(x=x,y=y[,i],family='gaussian',alpha=1,nlambda=50,standardize=FALSE,type.gaussian='covariance',penalty.factor=wlas)
			apl.bic <- log(colSums((predict(apl,x)-matrix(y[,i],nrow=nrow(y),ncol=ncol(apl$beta)))^2)/apl$nobs)
			apl.b	<- apl$beta[,which.min(apl.bic)]
			apl.a	<- apl$a0[which.min(apl.bic)]
			apl.vec	<-rep(0,ncol(aplas)-1)
			apl.vec[las.bmin[-1]!=0]	<-apl.b
			#print(c(sum(las.bmin[-1]!=0),sum(apl.b!=0),ncol(aplas)-1))
			
			#print(length(apl.b))
			#print(length(las.bmin[-1]!=0))
			aplas[i,]	<- c(apl.a,apl.vec)
		#	print(apl.b)
		#	print(las.b)
		#	print(ada.b)
		}
		else{aplas[i,1]<-mean(y[,i])}
	}
	adaest<-list('ada.b'=ada.b,'whichvar'=whichvar,'y'=y,'gamma'=gamma,'gnet'=ada.mod,'bols'=bols,'apl.b'=aplas,'las.b'=las.b)

	return (adaest)
}	

# Forecasting using SW2005 JBES approach:
pcaesti<-function(whichvar,sldata,nbrfac,nbrser,fchor,starty,nlag)
{
	pc	<-princomp(sldata)
	pca	<-pc$score[,c(1:nbrfac)]
	pca	<-ts(pca,start=starty,freq=12)
	bSW	<-NULL

	for(i in whichvar)
	{
		for(l in 1:nlag){x	<- ts.intersect(y[,i],lag(y[,i],1-l-fchor))}
		x	<-ts.intersect(x,lag(pca,1-fchor))

		bpc	<-lm.fit(y=x[,1],x=cbind(1,x[,-1]))$coefficients
		bSW	<-rbind(bSW,bpc)
	}
	SWest<-list('bSW'=bSW,'nbrser'=nbrser,'whichvar'=whichvar,'pca'=pca,'nlag'=nlag,'fchor'=fchor)

return(SWest)
}
# As before, but using BIC to select nbr factor and lag order
pcaestiBIC<-function(whichvar,sldata,maxfac,nbrser,fchor,starty,maxlag)
{
	pc	<-princomp(sldata)
	pca	<-pc$score[,c(1:maxfac)]
	pca	<-ts(pca,start=starty,freq=12)
	BICSW	<-matrix(0,ncol=maxfac,nrow=maxlag)
	bSW	<-NULL

	for(l in 1:maxlag)
	{
		for(p in 1:maxfac)
		{
			ppca	<-pca[,1:p]
			SSR	<-mclapply(whichvar,mcSWreg,sldata,ppca,l)
			BICSW[l,p]	<-log(Reduce('+',SSR))/length(sldata-fchor-l-p) + (l+p)*log(length(sldata-fchor-l-p))/length(sldata-fchor-l-p)
		}
	}
	min.ind	<-which(BICSW==min(BICSW),arr.ind=TRUE)
	#print(min.ind)
	for(i in whichvar){
		y	<-sldata[,i]
		for(m in 1:min.ind[1]){y<-ts.intersect(sldata[,i],lag(y,-m-fchor+1))}
		xBIC	<-ts.intersect(y,pca[,1:min.ind[2]])

		bpc	<-lm.fit(y=xBIC[,1],x=cbind(1,xBIC[,-1]))$coefficients
	#	print(bpc)
		bSW	<-rbind(bSW,bpc)
	}
#	print(BICSW)

	SWest<-list('bSW'=bSW,'nbrser'=nbrser,'whichvar'=whichvar,'pca'=pca[,c(1:min.ind[2])],'nlag'=min.ind[1],'BIC.cf'=min.ind[2],'fchor'=fchor)

return(SWest)
}

mcSWreg	<-function(i,sldata,ppca,l)
{
#	print(i)
	y	<-sldata[,i]
	for(m in 1:l){y<-ts.intersect(y,lag(sldata[,i],-m-fchor+1))}
	x	<-ts.intersect(y,ppca)
#	print(dim(x))
#	print(x[c(1,2),-1])
	r<-lm.fit(y=x[,1],x=cbind(1,x[,-1]))$residuals		
	return(sum(r^2))
}


pcafc<-function(y.fc,y.est,SWest)
{
	fc.fit	<-NULL
	fc.err	<-NULL

	for(i in 1:length(SWest$whichvar))
	{
	#	print(SWest$pca)
	#	print(rhs)
		rhs	<-NULL
		for(l in 1:SWest$nlag){rhs<- ts.intersect(rhs,lag(y.est[,SWest$whichvar[i]],1-l-fchor))}
		rhs	<-ts.intersect(rhs,lag(SWest$pca,1-SWest$fchor))

	#	print(c(1,rhs[nrow(rhs),]))
#		print(SWest$nlag)
#		print(SWest$BIC.cf)
#		print(SWest$bSW)
#		print(ncol(rhs))
		fit	<- c(1,rhs[nrow(rhs),])%*%matrix(SWest$bSW[i,],ncol=1)
		fc.fit	<- cbind(fc.fit,fit)
		err	<- y.fc[nrow(y.fc),SWest$whichvar[i]]-fit
		fc.err	<- cbind(fc.err,err)
	}


	#SANITY CHECK:
	se	<-sqrt(diag(var(y.est[,1:131])))	
	insane	<-(fc.fit>(y.fc[1:131]+3*se))||(fc.fit<(y.fc[1:131]-3*se))
	fc.fit[insane]	<-NA
	fc.err[insane]	<-NA


		#print(fc.fit)
		#print(y.fc[nrow(y.fc),SWest$whichvar])
		#print(fc.err)

		SWfc<-list('fit'=fc.fit,'err'=fc.err,'rhs'=rhs)
		return(SWfc)
}




adafc<-function(esti,y.fc,nbrser,whichvar,sldata,fchor)
{
fc.err	<-NULL
fc.fit	<-NULL
nlag	<-ncol(esti[,-1])/nbrser

y	<-sldata
if(nlag>1){for(l in 1:(nlag-1)){y<-ts.intersect(y,lag(sldata,-l))}}

#	if(nlag==1){
#		par.fc	<-cbind(esti[,1],(esti[,-1] %^% fchor)) #this might happen to be plain wrong...
#		vec.fc	<-c(1,y[nrow(y),])%*%t(matrix(par.fc,ncol=nbrser+1))
#		}
	if(nlag>0){
		y.last	<- y[nrow(y),] #init last
		for(t in 1:fchor){
			y.tmp	<- esti%*%c(1,y.last) #forecasting next period
			if(nlag>1)y.last	<- c(y.tmp,y.last[c(1:((nlag-1)*nbrser))])
			if(nlag==1)y.last	<-y.tmp
			}
		vec.fc	<-y.tmp
		}


	for(i in c(1:length(whichvar)))
	{
		fc	<- vec.fc[whichvar[i]]
		fc.fit	<- cbind(fc.fit,fc)
		fc.err	<- cbind(fc.err,y.fc[fchor,whichvar[i]]-fc)
	}
colnames(fc.err)	<-colnames(y.fc)[whichvar]
colnames(fc.fit)	<-colnames(y.fc)[whichvar]

	#SANITY CHECK:
	se	<-sqrt(diag(var(y[1:131])))	
	insane	<-(fc.fit>(y.fc[1:131]+3*se))||(fc.fit<(y.fc[1:131]-3*se))
	fc.fit[insane]	<-NA
	fc.err[insane]	<-NA


fc.err<-list('err'=fc.err,'fit'=fc.fit)
	    
return(fc.err)
}
