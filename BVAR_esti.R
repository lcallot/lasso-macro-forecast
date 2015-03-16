##
##	A function implementing the BVAR estimation by Banbura Giannone and Reichlin 2010 JAE
##
##

BGR_BVAR <- function(Y,X,sigma,delta,nlag,lambda)
{

	epsilon	 <- 0.1

	# Constructing the dummy observation to match the minnesota prior:

	sigdel	 <- diag(sigma*delta)
	Yp  	 <- rbind(sigdel/lambda,matrix(0,length(sigma)*(nlag-1),length(sigma)),diag(sigma),matrix(0,1,length(sigma)))
	Xp	 <-  cbind(rbind(matrix(0,ncol=1,nrow=length(sigma)*(nlag+1)),epsilon),rbind(diag(1:nlag)%x%diag(sigma)/lambda,matrix(0,length(sigma)+1,length(sigma)*nlag)))


#	print(Yp)
#	print(Xp)
#	print(dim(Y))
#	print(dim(X))

	Yf <- rbind(Y,Yp)
	Xf <- rbind(X,Xp)

#	print(dim(Yf))
#	print(dim(Xf))

	bgr <- ((lm(Yf~Xf-1)))

#	print(det(cov(bgr$res)))
	BIC <-	-2*as.double(determinant(var(bgr$res))$modulus) 
	bgr$bic<-BIC
#	print(lambda)
	bgr$lambda=lambda
#	browser()
	return(bgr)
}



bvar.est<-function(yd.s.est,nbrser)
	{
	ols.var.d	<-Reduce('rbind',mclapply(1:nbrser,function(i,y,x)
				{lm.fit(y=y[,i],x=x)->z;return(c(z$coef,var(z$res)))}
				,yd.s.est,cbind(1,yd.s.est[,-c(1:nbrser)])))

	# saving variances and parameters
	sigma		<- ols.var.d[,ncol(ols.var.d)]
	ols.var.d	<- ols.var.d[,-ncol(ols.var.d)]
	fc.olsd		<-adafc(ols.var.d,y.fc,nbrser,whichvar,yd.s.est[,c(1:nbrser)],1)



#	print(length(sigma))
	xlag	<-ncol(yd.s.est)/nbrser-1
	#sigma	<-rep(sigma,xlag)	
	delta 	<- rep(0,length(sigma))
#	print(xlag)

	# SEARCH
	bvar.0<-BGR_BVAR(yd.s.est[,1:nbrser],cbind(1,yd.s.est[,-c(1:nbrser)]),sigma,delta,nlag=xlag,lambda=0)
	bvar.1<-BGR_BVAR(yd.s.est[,1:nbrser],cbind(1,yd.s.est[,-c(1:nbrser)]),sigma,delta,nlag=xlag,lambda=1)
	lambda0	<-10^-5
	lambda1	<-10
	iter	<-1

	while(abs(bvar.0$bic-bvar.1$bic)>0.001){
		if(bvar.0$bic>bvar.1$bic){
					lambda0	<-(lambda1-lambda0)/2
					bvar.0	<-BGR_BVAR(yd.s.est[,1:nbrser],cbind(1,yd.s.est[,-c(1:nbrser)]),sigma,delta,nlag=xlag,lambda=lambda0)
				}	
				else{
					lambda1	<-lambda0+(lambda1-lambda0)/2
					bvar.1	<-BGR_BVAR(yd.s.est[,1:nbrser],cbind(1,yd.s.est[,-c(1:nbrser)]),sigma,delta,nlag=xlag,lambda=lambda1)

			}
			iter<-iter+1
		}
	print(c(bvar.0$bic,bvar.1$bic,lambda0,lambda1,iter-1))		
	bvar		<-bvar.0
	bvar$coef	<-t(bvar$coef)
	bvar$ols2	<-ols.var.d
	return(bvar)
	}
