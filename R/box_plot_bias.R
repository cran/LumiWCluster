box_plot_bias <-
function(x,y,xlabName,ylabName,mainName){
	plot(factor(x),y,xlab=xlabName,ylab=ylabName,main=mainName,cex.lab=1,cex.axis=1,cex.main=1.5)
	points(factor(lowess(x,y)$x),lowess(x,y)$y,col='gray',type='l',lwd=2)
}

