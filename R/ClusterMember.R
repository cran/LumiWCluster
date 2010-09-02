ClusterMember <-
function(LumiWCluster,Sample_ID){

	cat('Optimal number of clusters:', LumiWCluster$opt.K,'\n')

	for(k in 1:LumiWCluster$opt.K){
		cat('Sample IDs in Cluster',k,'\n')	
		cat(as.character(Sample_ID[which(LumiWCluster$clusterID==k)]),'\n')
	}
}

