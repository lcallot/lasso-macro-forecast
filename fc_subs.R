
library('reshape2')
library('ggplot2')
library('scales')


plot.rmseT<-function(rmseT)
	{
	pT<- ggplot(rmseT,aes(colour=variable,linetype=variable,y=value,x=Date)) + geom_line() 
	pT<- pT	+ facet_wrap(~variable, scales="free_y") + theme_bw()
#	pT<- pT	+ scale_x_date(breaks = "2 year", minor_breaks = "1 year", labels=date_format("%Y"))
	ggsave(pT,file='plots/rmseT.pdf')
	}


rmse.vec=function(est.err,data)
	{
	scl.dat	=std.rmse(est.err,data)
	return(rowMeans(scl.dat^2))
	}



std.rmse<-function(rmseT,data)
	{
	sig.vec=diag(var(data))
	sig.mat=matrix(rep(sqrt(sig.vec),nrow(rmseT)),ncol=ncol(rmseT),byrow=TRUE)

	scl.dat=rmseT/sig.mat
	return(scl.dat)
	}
	

plot.rmseT.grp<-function(rmseT)
	{
	merp<- ggplot(mer2,aes(colour=variable,linetype=variable,y=value,x=Date)) + geom_line() 
	merp<-merp+ theme_bw() 
#	merp<-merp+scale_x_date(breaks = "2 year", minor_breaks = "1 year", labels=date_format("%Y")) 
	ggsave(merp,file='plots/merp.pdf')
	}


# Write a small routine to store the variable and equation index of parameters that are in the xth quantile of the of the parameters distributions.

big.param<-function(param,quant)
{
	x	<-melt(data.frame(abs(param)))
	cut	<-quantile(x$value,quant)
	print(paste('The cut of point is ',cut))
	big.ind	<-which(param>cut,arr.ind=TRUE)

return(big.ind)
}















#	1:	for a given param mat: return properly formated freq data

freq.quan<-function(param,quant,par.name,nbr.iter)
	{

	fq	<-list()
	nbr.eq	<-nrow(param)/nbr.iter
	

	param	<-abs(param[,-1]) #removing the constant
	par.nz	<-param[param!=0]

	cut	<-quantile(par.nz,quant)
	cut	<-0
	fq$cut	<-cut
	big.ind	<-NULL
	for(l in 1:2){for(i in 1:nbr.iter){big.ind	<-rbind(big.ind,which(param[(1+(i-1)*nbr.eq):(i*nbr.eq),(1+(l-1)*nbr.eq):(l*nbr.eq)]>cut,arr.ind=TRUE))}}
	colnames(big.ind)<-c('Equation','Covariate')	
	bog.ind	<-big.ind[order(big.ind[,1],big.ind[,2]),]

	#construct a matrix of eq/cov frequency
	freq.mat	<-table(factor(big.ind[,1],1:nbr.eq),factor(big.ind[,2],1:ncol(param)/2))
	long.frq	<-melt(freq.mat)
	names(long.frq)	<-c(colnames(big.ind),'count')

	fq$frq	<-long.frq
	fq$iter	<-nbr.iter
	#browser()
	return(fq)
	}



# Averages the parameters across forecast horizon
size.par<-function(param,par.name,nbr.iter)
	{

	fq	<-list()
	nbr.eq	<-nrow(param)/nbr.iter
	

	param	<-abs(param[,-1]) #removing the constant

	fq$cut	<-0
	big.ind	<-matrix(0,nrow=nbr.eq,ncol=ncol(param))
	for(i in 1:nbr.iter){big.ind	<-big.ind+param[(1+(i-1)*nbr.eq):(i*nbr.eq),]}
	big.ind	<-big.ind/nbr.iter
	rownames(big.ind)<-1:nbr.eq
	colnames(big.ind)<-1:ncol(param)

	#construct a matrix of eq/cov frequency
	long.frq	<-melt(big.ind)
	names(long.frq)	<-c('Equation','Covariate','count')

	fq$frq	<-long.frq
	fq$iter	<-nbr.iter
	#browser()
	return(fq)
	}

#	2:	agg the data for all esti, mk heat map
agg.freq<-function(quant,par.agg,nbr.iter)
	{
	agg	<-list()
	agg$par	<-NULL
	agg$cut	<-NULL
	for(p in names(par.agg))
	{
		#fq	<-freq.quan(par.agg[[p]],quant,p,nbr.iter)
		fq	<-size.par(par.agg[[p]],p,nbr.iter)
		agg$par	<-rbind(agg$par,cbind(rep(p,nrow(fq$frq)),fq$frq))
		agg$cut	<-c(agg$cut,fq$cut)
	}
	colnames(agg$par)[1]<-'Estimator'
	agg$quant<-quant
	agg$iter <-fq$iter
	return(agg)
	}

plot.agg.freq<-function(agg.count)
{
	levels(agg.count$par[[1]])<-paste(levels(agg.count$par[[1]]),'\n Cut-point: ',round(agg.count$cut,1),sep='')	
	x_lab<-seq(1,131,10)
	y_lab<-seq(1,262,10)
	agg.count$par$count[agg.count$par$count==0]<-NA

 	frq.plot<- ggplot(agg.count$par, aes(Equation,Covariate)) + geom_tile(aes(fill = count),colour =   "white") 
	frq.plot<- frq.plot + facet_wrap(~Estimator)
	frq.plot<- frq.plot + scale_x_discrete(breaks=x_lab)
	frq.plot<- frq.plot + scale_y_discrete(breaks=y_lab)
 	frq.plot<- frq.plot + scale_fill_gradient2(na.value='white',low='white', high = "navyblue")+theme_bw()
#	frq.plot<- frq.plot + ggtitle(paste('Count of parameters above the ',100*agg.count$quant,'% quantile (of non-zero parameters)\n',agg.count$iter,' forecasts',sep=''))
#	print(frq.plot)
	ggsave(frq.plot,file=paste('plots/avg_param','','.pdf',sep=''))
}


plot.agg.indep<-function(agg.count)
{
	x_lab<-seq(1,131,10)
	y_lab<-seq(1,262,10)
#	agg.count$par$count[agg.count$par$count==0]<-NA

	pdf('plots/ind_avg_param.pdf',width=14,height=21)
	i=1
	for(est in levels(agg.count$par$Estimator))
		{
		par.est<-agg.count$par[agg.count$par$Estimator==est,]
		## set color representation for specific values of the data distribution
		quantile_range 	<- quantile(agg.count$par$count[agg.count$par$Estimator==est], probs = seq(0, 1, 0.2),na.rm=TRUE)
		## prepare label text (use two adjacent values for range text)
		label_text 	<- rollapply(round(quantile_range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
		## use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
		color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range) - 1)

		frq.plot<- ggplot(data=subset(agg.count$par,Estimator==est),aes(Equation,Covariate)) + geom_tile(aes(fill = count),colour =   "white") +theme_bw()
		frq.plot<- frq.plot + scale_x_discrete(breaks=x_lab)
		frq.plot<- frq.plot + scale_y_discrete(breaks=y_lab)
		#frq.plot<- frq.plot + scale_fill_gradient2(na.value='white',low='red', high = "navyblue")
		frq.plot<- frq.plot + scale_colour_manual(values = color_palette, name = "",breaks=quantile_range,labels = label_text) 
		frq.plot<- frq.plot + ggtitle(est)
		assign(paste('p',i,sep=''),frq.plot)
		i<-i+1
		}

	ind.plot<-arrange(p1,p2,p3,p4,p5,p6,p7,ncol=2)
#	ind.plot<- ind.plot + ggtitle(paste('Count of parameters above the ',100*agg.count$quant,'% quantile (of non-zero parameters)\n',agg.count$iter,' forecasts',sep=''))

	print(ind.plot)
	dev.off()
#	ggsave(ind.plot,file=paste('plots/ind_avg_param','','.pdf',sep=''))
}




#Aggregate the location of the params greater than a given quantile for each estimator, vizualization...
big.loc<-function(adaest,quant)
	{
	par.names<-c('las.b','ada.b','apl.b')
	all.cut	<-NULL
	all.loc	<-NULL
	for(p in par.names)
	{
		#cut	<-quantile(adaest[[p]],quant)
		cut	<-100
		big.ind	<-which(adaest[[p]]>cut,arr.ind=TRUE)
		all.cut	<-c(all.cut,round(cut,1))
		rownames(big.ind)<-rep(p,nrow(big.ind))
		all.loc	<-rbind(all.loc,big.ind)
	}
	colnames(all.loc)<-c('equation','covariate')
	all.cut<-matrix(all.cut,ncol=1)
	rownames(all.cut)<-par.names
	colnames(all.cut)<-'quantile'
	
	#construct a matrix of eq/cov frequency
	freq.mat<-table(all.loc[,1],all.loc[,2])
	long.frq<-melt(freq.mat)
	names(long.frq)<-c(colnames(all.loc),'count')
	#data.m	<- ddply(long.fr, .(), transform, rescale = rescale(value))

	
 	frq.plot<- ggplot(long.frq, aes(equation,covariate)) + geom_tile(aes(fill = count),colour =   "white") 
 	frq.plot<- frq.plot + scale_fill_gradient(low = "white", high = "steelblue")+theme_bw()
	#frq.plot<- frq.plot + annotate("table", x=132, y=200, table=all.cut, just=c("left", "top"))#,
             #theme=theme.list(show.box = TRUE, separator = "black",
              #show.csep = TRUE, show.rsep = TRUE, show.colnames=T))

#create the grob:
	tab.grob<-tableGrob(all.cut,theme=theme_bw(),show.box=TRUE,gpar.coretext = gpar(fontsize=8,col="black"),gpar.coltext = gpar(fontsize=10,col="black",fontface='bold'),gpar.rowtext = gpar(fontsize=8,col="black"), gpar.colfill = gpar(fill=NA,col=NA), gpar.rowfill = gpar(fill=NA,col=NA),gpar.corefill=gpar(fil=NA,col=NA))

	frq.plot<-frq.plot+annotation_custom(grob=tab.grob,xmin=137,ymax=50)

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(frq.plot))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#	print(frq.plot)
	}
















#return the matrix of quantiles of the lasso estimators for a run 
quant.par<-function(adaest,quant)
{
	par.names<-c('las.b','ada.b','apl.b')
	all.q	<-NULL
	for(p in par.names)
	{
		par.q	<-quantile(matrix(adaest[[p]],ncol=1),quant)
		all.q	<-cbind(all.q,par.q)
	}
	colnames(all.q)<-par.names
return(all.q)
}

#plot the density of the params of the lasso estimators for a run 
dens.par<-function(adaest)
{	
	par.names<-c('las.b','ada.b','apl.b')
	all.par	<-NULL
	for(p in par.names)
	{
		all.par<-cbind(all.par,matrix(adaest[[p]],ncol=1))
	}
	colnames(all.par)<-par.names
	long.par	<-melt(data.frame(all.par))
	

	# DENSITY UNREADABLE, long tails...
	#dens.par	<-ggplot(long.par,aes(x=value,colour=variable))+geom_density()
	#print(dens.par)

}
