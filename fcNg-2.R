library("dse")
library('expm')
library("lars")
library('glmnet')
library("multicore")
library("xtable")
library('ggplot2')
library('reshape2')

source("SerenaData.R")
source("ShrinkNg.R")
source('BVAR_esti.R')
source('fc_subs.R')

sldata	<-ts(read.table('SerenaDat'),start=c(1964,3),freq=12)



for(fchor in c(3,6,12)){ #horizon selection
	starty	<-	c(1964,3)
	endy	<-	c(1999,12)
	fsty	<-	incdate(c(1999,12),fchor)
	maxt	<-	12*(2007-endy[1])+(12-endy[2])-fchor+1
	print(paste('nbr fc computed: ',maxt,sep=''))

	maxlag	<-	2 #Setting max nbr lags
	gamma	<-	c(1)

	whichvar<-c(1:131)

	# GROUPING ALONG THE LINE OF LUDVIGSON AND NG
	blanktype			<-rep(0,131)
	blanktype[c(1,2,6:20)]		<-1
	blanktype[c(21:50,129:131)]	<-2
	blanktype[c(51:60)]		<-3
	blanktype[c(61:70,3:5,132)]	<-4
	blanktype[c(71:81)]		<-5
	blanktype[c(86:107)]		<-6
	blanktype[c(108:128)]		<-7
	blanktype[c(82:85)]		<-8
	#for some reason series 46 isn't in the data set.
	LNtype	<-blanktype[-46]
	LNnames	<-c('Out and Inc','Labor mkt','Housing','Cons, Ord, Inv','Money and cred','Bond and XR','Prices','Stock mkt')

	est.names	<-c('SWB','SW1','SW3','SW5','olsd','adad','lasd','grpd','bvar','apla','ols2')
	na.names	<-paste(est.names,'.na',sep='')
	err.names	<-paste(est.names,'.err',sep='')

	for(n in err.names){assign(n,NULL)}


	sldata	<-	window(sldata,start=starty,end=c(2007,12))
	nbrser	<-	ncol(sldata)
	yd	<-	sldata
	yd.s	<-	sldata
	y	<-	sldata
	y.s	<-	sldata

	print('data loaded, crunching time')
	#Lag selection for VAR
	#lag.dir <-VARBIC(sldata,maxlag,direct=TRUE,fchor)
	lag.dir<-1

	#For the ols lag
	for(l in 1:lag.dir){yd	<-ts.intersect(yd,lag(sldata,-l-fchor+1))}
	for(l in 1:lag.dir){y	<-ts.intersect(y,lag(sldata,-l-fchor+1))}

	#for the Lasso and bayesian VAR 
	for(l in 1:maxlag){yd.s	<-ts.intersect(yd.s,lag(sldata,-l-fchor+1))}
	for(l in 1:maxlag){y.s	<-ts.intersect(y.s,lag(sldata,-l-fchor+1))}

	obs	<-nrow(y.s)	
	if(!(nbrser*maxlag<obs)) stop(paste("Less observations than variables, I don't like that\n obs:",obs,'param:',nbrser*maxlag,sep=''))
	print(paste('obs:',obs,sep=''))
	print(paste('param per eq: ',nbrser*maxlag+1,sep='')) 


	selfreq	<-list('ada.dir'=NULL,'las.dir'=NULL,'grp.dir'=NULL,'aplas'=NULL)
	par.agg <-list('Lasso'=NULL,'aLasso (ols)'=NULL,'aLasso (lasso)'=NULL,'Group'=NULL,'VAR(1)-OLS'=NULL,'BVAR'=NULL,'VAR(2)-OLS'=NULL)
#maxt=9
	for(t in c(1:maxt))
		{
		#House keeping
		print(paste('Sample',t,sep=''))
		y.est	<-window(y,end=endy)
		yd.est	<-window(yd,end=endy)
		y.fc	<-window(y,start=incdate(endy,1),end=incdate(endy,fchor))
		nbrser	<-length(colnames(sldata))

		y.s.est	<-window(y.s,end=endy)
		yd.s.est<-window(yd.s,end=endy)
		y.s.fc	<-window(y.s,start=incdate(endy,1),end=incdate(endy,fchor))


		print('VAR')
		y.var		<-yd.est
		#Estimating the VAR by OLS
		ols.var.d	<-Reduce('rbind',mclapply(1:nbrser,function(i,y,x){lm.fit(y=y[,i],x=x)->z;return(c(z$coef,var(z$res)))},y.var,cbind(1,y.var[,-c(1:nbrser)])))

		# saving variances and parameters
		sigma		<- ols.var.d[,ncol(ols.var.d)]
		ols.var.d	<- ols.var.d[,-ncol(ols.var.d)]
		fc.olsd		<-adafc(ols.var.d,y.fc,nbrser,whichvar,y.var[,c(1:nbrser)],1)

		print('Lasso and adaptive Lasso')
		adaestd		<-mcadaesti(whichvar,yd.s.est,gamma,nbrser)
		selfreq$las.dir	<-cbind(selfreq$las.dir,adaestd$las.b[,-1]!=0)
		selfreq$ada.dir	<-cbind(selfreq$ada.dir,adaestd$ada.b[,-1]!=0)
		selfreq$aplas	<-cbind(selfreq$aplas,adaestd$apl.b[,-1]!=0)
		selfreq$ols2	<-cbind(selfreq$ols2,adaestd$bols[,-1]!=0)


		print('BVAR')
		#delta <- rep(0,length(sigma))
		#bvar.d<-BGR_BVAR(y.var,cbind(1,y.var[,-c(1:nbrser)]),sigma,delta,nlag=lag.dir,lambda=0)
		#for(lambda in c(1:10)/10){
		#	bvar <- BGR_BVAR(y.var,cbind(1,y.var[,-c(1:nbrser)]),sigma,delta,nlag=lag.dir,lambda)
		#	if(bvar$bic<bvar.d$bic)bvar.d<-bvar#bic selection of lambda
		#}
		#if(length(bvar.d)==1)bvar.d=bvar  #Using the last BVAR if none of them was selected before.

		#bvar	<-bvar.est(yd.s.est,nbrser	
		
		bvar	<-bvar.est(yd.s.est,nbrser)
		fc.bvar	<-adafc(bvar$coef,y.fc,nbrser,whichvar,yd.s.est[,c(1:nbrser)],fchor)
		fc.ols2	<-adafc(bvar$ols2,y.fc,nbrser,whichvar,yd.s.est[,c(1:nbrser)],fchor)


		#Stock and Watson forecasts
		print('Common factor')
		SWest.names	<- paste('SWest',c(c(1,3,5),'BIC'),sep='')
		SWfc.names	<- paste('SWfc',c(c(1,3,5),'BIC'),sep='')

		for(cf in c(1,3,5)){
			assign(SWest.names[cf],pcaesti(whichvar,window(sldata,start=starty,end=endy),nbrfac=cf,nbrser,fchor,starty,nlag=1))
			assign(SWfc.names[cf],pcafc(y.fc,y.est,get(SWest.names[cf])))
			assign(paste('SW',cf,'.err',sep=''),rbind(get(paste('SW',cf,'.err',sep='')),get(SWfc.names[cf])$err))
		}
#		#SW esti using BIC
		SWestBIC<-pcaestiBIC(whichvar,window(sldata,start=starty,end=endy),maxfac=5,nbrser,fchor,starty,maxlag=1)
		SWfcBIC	<-pcafc(y.fc,y.est,SWestBIC)
		SWB.err <-rbind(SWB.err,SWfcBIC$err)

		# Group Lasso
		print('Group Lasso estimation')
	#	grp.vec	<-rep(LNtype,maxlag)
		grp.vec <-c(LNtype,LNtype+8)
		grpd	<-estgrp(grp.vec,nbrser,yd.s.est,gamma[1])
		selfreq$grp.dir	<-cbind(selfreq$grp.dir,grpd$grp[,-1]!=0)


		print('Forecasting with regularized estiamtors')
		fc.gBICd	<-adafc(adaestd$ada.b,y.s.fc,nbrser,whichvar,adaestd$y[,c(1:nbrser)],1)
		fc.lBICd	<-adafc(adaestd$las.b,y.s.fc,nbrser,whichvar,adaestd$y[,c(1:nbrser)],1)
		fc.apla		<-adafc(adaestd$apl.b,y.s.fc,nbrser,whichvar,adaestd$y[,c(1:nbrser)],1)
		fc.grpd		<-adafc(grpd$grp,y.s.fc,nbrser,whichvar,grpd$y[,c(1:nbrser)],1)
		


		#param aggregation
		par.agg[[1]]<-rbind(par.agg[[1]],adaestd$las.b)
		par.agg[[2]]<-rbind(par.agg[[2]],adaestd$ada.b)
		par.agg[[3]]<-rbind(par.agg[[3]],adaestd$apl.b)
		par.agg[[4]]<-rbind(par.agg[[4]],grpd$grp)
		par.agg[[5]]<-rbind(par.agg[[5]],ols.var.d)
		par.agg[[6]]<-rbind(par.agg[[6]],t(bvar$coeff))
		par.agg[[7]]<-rbind(par.agg[[7]],bvar$ols2)
		


		#Forecast error aggregation
		bvar.err	<-rbind(bvar.err,unclass(fc.bvar$err))
		ols2.err	<-rbind(ols2.err,unclass(fc.ols2$err))
		adad.err	<-rbind(adad.err,unclass(fc.gBICd$err))
		lasd.err	<-rbind(lasd.err,unclass(fc.lBICd$err))
		apla.err	<-rbind(apla.err,unclass(fc.apla$err))
		olsd.err	<-rbind(olsd.err,unclass(fc.olsd$err))
		grpd.err	<-rbind(grpd.err,unclass(fc.grpd$err))

	 #   starty  <-incdate(starty,1)
	    endy    <-incdate(endy,1)

	} # End of main loop!
	rownames(par.agg[4])<-NULL

	SWB.na		<-colSums(matrix(is.na(SWB.err),nrow=nrow(SWB.err)))
	SW1.na		<-colSums(matrix(is.na(SW1.err),nrow=nrow(SW1.err)))
	SW3.na		<-colSums(matrix(is.na(SW3.err),nrow=nrow(SW3.err)))
	SW5.na		<-colSums(matrix(is.na(SW5.err),nrow=nrow(SW5.err)))

	olsd.na		<-colSums(matrix(is.na(olsd.err),nrow=nrow(olsd.err)))
	ols2.na		<-colSums(matrix(is.na(ols2.err),nrow=nrow(ols2.err)))
	bvar.na		<-colSums(matrix(is.na(bvar.err),nrow=nrow(bvar.err)))
	adad.na		<-colSums(matrix(is.na(adad.err),nrow=nrow(adad.err)))
	lasd.na		<-colSums(matrix(is.na(lasd.err),nrow=nrow(lasd.err)))
	grpd.na		<-colSums(matrix(is.na(grpd.err),nrow=nrow(grpd.err)))
	apla.na		<-colSums(matrix(is.na(apla.err),nrow=nrow(apla.err)))

# Replacing insane forecasts with y_t+h-y_t

	last.obs	<-window(sldata,start=c(1999,12))
	last.err	<-diff(last.obs,fchor)  
	
	for(n in na.names)
		{
		x	<-get(n)
		x[is.na(x)]<-last.err[is.na(x)]
		}

	#SANE MSE	
	rmse.SWBIC	<-colMeans(SWB.err^2,na.rm=TRUE)
	rmse.SW1	<-colMeans(SW1.err^2,na.rm=TRUE)
	rmse.SW3	<-colMeans(SW3.err^2,na.rm=TRUE)
	rmse.SW5	<-colMeans(SW5.err^2,na.rm=TRUE)

	rmse.olsd    	<-colMeans(olsd.err^2,na.rm=TRUE)
	rmse.ols2    	<-colMeans(ols2.err^2,na.rm=TRUE)
	rmse.bvar	<-colMeans(bvar.err^2,na.rm=TRUE)
	rmse.adaBgd  	<-colMeans(adad.err^2,na.rm=TRUE)
	rmse.lasBgd  	<-colMeans(lasd.err^2,na.rm=TRUE)
	rmse.grpd    	<-colMeans(grpd.err^2,na.rm=TRUE)
	rmse.apla    	<-colMeans(apla.err^2,na.rm=TRUE)

	rmse		<-rbind(rmse.olsd,rmse.bvar,rmse.lasBgd,rmse.adaBgd,rmse.apla,rmse.grpd,rmse.ols2)
	rownames(rmse)	<-c('OLS-VAR','BVAR','Lasso','aLasso-OLS','aLasso-Lasso','agLasso','VAR(2)')

	rrmse	<-	rmse/matrix(rep(rmse.olsd,nrow(rmse)),ncol=ncol(rmse),byrow=TRUE)
	mse.SW	<-	rbind(rmse.SW1,rmse.SW3,rmse.SW5,rmse.SWBIC)
	rrmseSW	<-	mse.SW/matrix(rep(rmse.olsd,nrow(mse.SW)),ncol=ncol(rmse),byrow=TRUE)
	rownames(rrmseSW)	<-c('CF 1','CF 3','CF 5','CF BIC')
	#aggrr	<-	cbind(rowMeans(rrmse[,1:20]),rowMeans(rrmse[,21:50]),rowMeans(rrmse[,51:70]),rowMeans(rrmse[,71:81]),rowMeans(rrmse[,82:102]),rowMeans(rrmse[,103:107]),rowMeans(rrmse[,108:131]))


	aggrr	<-NULL
	aggrrSW	<-NULL
	aggVAR	<-NULL

	gmse.SW	<-NULL
	gmse.VAR<-NULL


	selgp	<-list('ada.dir'=matrix(NA,nrow=8,ncol=8),'las.dir'=matrix(NA,nrow=8,ncol=8),'apl.dir'=matrix(NA,nrow=8,ncol=8),'grp.dir'=matrix(NA,nrow=8,ncol=8))
	nbrgp	<-list('ada.dir'=matrix(NA,nrow=8,ncol=8),'las.dir'=matrix(NA,nrow=8,ncol=8),'apl.dir'=matrix(NA,nrow=8,ncol=8),'grp.dir'=matrix(NA,nrow=8,ncol=8))

	for(i in 1:8){
		aggrr	<-cbind(aggrr,rowMeans(rrmse[,LNtype==i]))
		aggrrSW	<-cbind(aggrrSW,rowMeans(rrmseSW[,LNtype==i]))
		aggVAR	<-cbind(aggVAR,mean(rmse.olsd[LNtype==i]))

		gmse.SW	<-cbind(gmse.SW,rowMeans(mse.SW[,LNtype==i]))
		gmse.VAR<-cbind(gmse.VAR,rowMeans(rmse[,LNtype==i]))
		for(j in 1:8){
			for(e in 1:length(selgp)){
				selgp[[e]][i,j]<-mean(colMeans(selfreq[[e]][LNtype==i,])[LNtype==j])
				nbrgp[[e]][i,j]<-sum(colSums(selfreq[[e]][LNtype==i,])[LNtype==j])
			}
		}
	}

	colnames(gmse.SW)		<-LNnames	
	colnames(gmse.VAR)		<-LNnames	
	rownames(gmse.SW)		<-rownames(rrmseSW)
	rownames(gmse.VAR)		<-rownames(rmse)

	print(paste('forecasting horizon ',fchor,sep=''))
	print('BIC selected lag order:')
	print(paste('nbr fc computed: ',maxt,sep=''))
	colnames(aggrr)		<-LNnames	
	print(aggrr)
	colnames(aggrrSW)	<-LNnames	
	print(aggrrSW)
	colnames(aggVAR)	<-LNnames	

###########################################
#		Should output be saved?
###########################################	

	if(TRUE){
		print(xtable(aggrr,digits=3,caption=paste('fc horizon:',fchor,' - nbr fc:',maxt,sep='')),file=paste('tabs/RR_fc_h_',fchor,sep=''),only.content=TRUE,append=FALSE,booktabs=TRUE)
		print(xtable(aggrrSW,digits=3,caption=paste('BIC lag: =',lag.dir,sep=' ')),file=paste('tabs/RR_fc_h_',fchor,sep=''),only.content=TRUE,append=TRUE,booktabs=TRUE)

	#	print(xtable(aggVAR,digits=4,caption='VAR MSE'),file=paste('new_fc_h_',fchor,sep=''),only.content=TRUE,append=TRUE,booktabs=TRUE)

	#	write('MSE all models (not relative)',file=paste('new_fc_h_',fchor,sep=''),append=TRUE)
	#	write(gmse.VAR,file=paste('new_fc_h_',fchor,sep=''),append=TRUE)
	#	write(gmse.SW,file=paste('new_fc_h_',fchor,sep=''),append=TRUE)



#WRITING THE SELECTION RATES!!!	
		write('Group selection RATE',file=paste('tabs/sel_fc_h_',fchor,sep=''),append=FALSE)
		for(e in 1:length(selgp)){print(xtable(selgp[[e]],digits=3,caption=names(selgp)[e]),file=paste('tabs/sel_fc_h_',fchor,sep=''),append=TRUE,only.content=FALSE)}

#NEED TO SCALE BY NUMBER OF EQUATIONS* NUMBER OF FORECASTS
		write('Group selection COUNT',file=paste('tabs/nbr_fc_h_',fchor,sep=''),append=FALSE)
		for(e in 1:length(nbrgp)){print(xtable(nbrgp[[e]],digits=3,caption=names(nbrgp)[e]),file=paste('tabs/nbr_fc_h_',fchor,sep=''),append=TRUE,only.content=FALSE)}



# SAVING AGG INSANCE FC	

		write('INSANE FORECASTS',file=paste('tabs/insane_h_multi',fchor,sep=''),append=FALSE)
		for(n in na.names)
		{
			x	<-get(n)
			xx	<-NULL
			for(i in 1:8){xx<-c(xx,mean(x[LNtype==i]))}
			write(n,append=TRUE,file=paste('tabs/insane_h_multi',fchor,sep=''))
			write(xx,append=TRUE,file=paste('tabs/insane_h_multi',fchor,sep=''))
		}
}



#Large parameters:
	agg.par<-agg.freq(0.99,par.agg,maxt)
	plot.agg.freq(agg.par)



#RMSE across time	
	relmseT<-data.frame(ts(cbind(rowMeans((adad.err/olsd.err)^2),rowMeans((apla.err/olsd.err)^2),rowMeans((lasd.err/olsd.err)^2),rowMeans((grpd.err/olsd.err)^2),rowMeans((SW1.err/olsd.err)^2)),start=fsty,end=incdate(fsty,maxt),freq=12))


	rmseT=NULL
	for(e in c('olsd','adad','lasd','grpd','SW1','bvar','ols2')){	rmseT<-cbind(rmseT,rmse.vec(get(paste(e,'.err',sep='')),sldata));print(dim(rmseT))}
	
	erm	<-cbind(seq.Date(as.Date('2000-01-01'),by='1 month',length.out=nrow(rmseT)),rmseT)
	colnames(erm)<-c('Date','VAR-OLS','aLasso.OLS','Lasso','agLasso','SW1','BVAR','VAR(2)')
	rownames(erm)<-NULL

	vrrr	<-melt(data.frame(erm),id=c(1))	
	plot.rmseT(vrrr)

	mer2<-melt(data.frame(erm),id=1)
	plot.rmseT.grp(mer2)


}	    



