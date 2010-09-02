LumiWCluster <-
function(logitbeta,w,g,K=c(2:5),lambda=exp(log(seq(1,50,length=20))),adaptive=TRUE,center=TRUE,normalize=FALSE,epsilon=1e-4,max.iter=2000,trace=TRUE){
	
	
	if (any(!is.finite(logitbeta))) 
		stop("infinite or missing values in 'logitbeta'")
	logitbeta <- t(logitbeta)
	n <- dim(logitbeta)[1]
	p <- dim(logitbeta)[2]
	ori.logitbeta <- logitbeta
	
	K <- unique(as.integer(K))  ## K no. of clusters -- positive integer greater than 1
	K <- K[K>0]
	if(length(K)>0){
		if(min(K)==1) K <- K[-which(K==1)]
	}
	if (missing(K)|length(K)==0) {
		K <- c(2:5)
	}
	# unique.logitbeta = unique(logitbeta)
	# nn <- nrow(unique.logitbeta)
	# if(any(K>nn)) {
	#	stop("There are ", nn, " unique observations. 'K' may not exceed this number.")
	# }
	cat('Searching for K=',K,'...','\n')
	maxK <- max(K)
	nK <- length(K)
	
	lambda <- unique(lambda[lambda>0])  ## lambda, penalty parameter, greater than 0
	if (missing(lambda)|length(lambda)==0) {
		lambda <- exp(log(seq(1,100,length=50)))
	}
	nLambda <- length(lambda)

	if (any(w<0)) 
		stop("Negative weights in 'w'")
	if (missing(w)){
		w <- rep(1,n)
	} 
	w <- w/sum(w)
	if(length(w)!=n) stop('Not all observations assigned weights')

	if (any(g < 0)) 
		stop("Negative weights in 'g'")
	if (missing(g)){
		g <- rep(1,p)
	} 
	g <- g/sum(g)
	if(length(g)!=p) stop('Not all CpG/genes assigned weights')

	
	mean.logitbeta <- rep(0,p)
	norm.logitbeta <- rep(1,p)
	
	if (center) {
		mean.logitbeta <- apply(logitbeta,2,mean)
		logitbeta <- scale(logitbeta,mean.logitbeta,FALSE)
	}

	if (normalize) {
		norm.logitbeta <- apply(logitbeta,2,sd)
		logitbeta <- scale(logitbeta,FALSE,norm.logitbeta)
	}

	zz <- .C("LumiWCluster", as.double(as.vector(t(ori.logitbeta))), as.double(as.vector(t(logitbeta))), as.integer(n), as.integer(p), as.double(mean.logitbeta), as.double(norm.logitbeta), as.double(w), as.double(g), as.integer(K), as.integer(nK), as.integer(maxK), as.double(lambda), as.integer(nLambda), as.integer(adaptive), as.double(epsilon), as.integer(max.iter), as.integer(trace), BIC=double(nK*nLambda), opt.BIC=double(1), opt.K=integer(1), opt.lambda=double(1), ClusterID=integer(n), converge = integer(1), pi.vec=double(maxK), mu=double(maxK*p), mu.centering=double(maxK*p), init.mu = double(maxK*p),  sigma2.vec=double(p), iter=integer(1), package="LumiWCluster")

	opt.K <- zz$opt.K
	opt.lambda <- zz$opt.lambda
	BIC <- zz$BIC
	BIC <- matrix(BIC, ncol=nLambda, byrow=T)
	opt.BIC <- zz$opt.BIC
	clusterID <- zz$ClusterID
	converge <- zz$converge
	pi.vec <- zz$pi.vec
	pi.vec <- pi.vec[1:opt.K]
	mu <- zz$mu
	mu <- matrix(mu, ncol=p, byrow=T)
	mu <- mu[1:opt.K,]
	mu.centering <- zz$mu.centering
	mu.centering <- matrix(mu.centering, ncol=p, byrow=T)
	mu.centering <- mu.centering[1:opt.K,]
	init.mu <- zz$init.mu
	init.mu <- matrix(init.mu, ncol=p, byrow=T)
	init.mu <- init.mu[1:opt.K,]
	sigma2.vec <- zz$sigma2.vec
	iter <- zz$iter
		
	cat('DONE! Optimal K =',opt.K,'\n')
	
	return(list(opt.K=opt.K,opt.lambda=opt.lambda,clusterID=clusterID,BIC=opt.BIC,BIC_iter=BIC,lambda=lambda,K=K,converge=converge,pi.vec=pi.vec,mu=mu,mu.centering=mu.centering,init.mu=init.mu,sigma2.vec=sigma2.vec,iter=iter))

}

