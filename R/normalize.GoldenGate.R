normalize.GoldenGate <-
function(Probe_ID, beta, det.p.value, design.file, plot.option = FALSE){

	if(is.vector(beta) == TRUE) beta <- matrix(beta)
	if(is.vector(det.p.value) == TRUE) det.p.value <- matrix(det.p.value)

	beta[beta == 0] <- min(beta[ beta > 0],1e-3)/10
	beta[beta == 1] <- max(beta[ beta < 1],0.999) + (1 - max(beta[ beta < 1],0.999))/100

	matchID <- rep(NA,length(Probe_ID))

	for(i in 1:length(matchID)){
		matchID[i] <- which(as.character(design.file$Probe_ID) == as.character(Probe_ID[i]))
	}

	L <- design.file$L[matchID] 
	GC <- design.file$GC[matchID]

	y_raw <- log(t(beta)/(1 - t(beta)))

	n <- dim(y_raw)[1]
	p <- dim(y_raw)[2]
	y <- matrix(0,nrow=n,ncol=p)
	
	w <- apply(log(t(det.p.value)),1,median)
	w <- w/sum(w)
	g <- apply(log(t(det.p.value)),2,median)
	g <- g/sum(g)

	L_new <- L
	L_new[which(L<=44)] <- 44
	L_new[which(L>=57)] <- 57
	GC_new <- GC
	GC_new[which(GC<=0.4)] <- 0.4
	GC_new[which(GC>=0.8)] <- 0.8

	L_coef <- 0.04458
	GC_coef <- -3.658

	for(i in 1:n){
		
		y[i,] <- y_raw[i,] - L_coef*L_new - GC_coef*GC_new + mean(L_coef*L_new) + mean(GC_coef*GC_new)

	}
	
	norm.beta <- t(inv_logit(y))

	if(plot.option == TRUE ) plotBAnorm(beta, norm.beta, design.file)

	return(list(norm.beta = norm.beta, norm.logitbeta = t(y), w = w, g = g))

}

