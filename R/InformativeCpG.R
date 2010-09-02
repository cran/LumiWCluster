InformativeCpG <-
function(LumiWCluster,Probe_ID){

	nonZero <- which(apply(abs(LumiWCluster$mu.centering),2,sum) != 0)
	return(as.character(Probe_ID[nonZero]))
	
}

