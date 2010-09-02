inv_logit <-
function(mu){
	return(exp(mu)/(1+exp(mu)))
}

