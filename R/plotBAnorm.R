plotBAnorm <-
function(ori.beta, norm.beta, design.file){

	par(mfrow = c(2,2))
	box_plot_bias(rep(design.file$L,dim(ori.beta)[2]), c(ori.beta),xlabName='Sequence length',ylabName='Unnormalized beta',mainName='Before Norm')
	box_plot_bias(rep(round(design.file$GC*20)/20,dim(ori.beta)[2]), c(ori.beta),xlabName='GC content',ylabName='Unnormalized beta',mainName='Before Norm')
	box_plot_bias(rep(design.file$L,dim(ori.beta)[2]), c(norm.beta),xlabName='Sequence length',ylabName='Normalized beta',mainName='After Norm')
	box_plot_bias(rep(round(design.file$GC*20)/20,dim(ori.beta)[2]), c(norm.beta),xlabName='GC content',ylabName='Normalized beta',mainName='After Norm')


}

