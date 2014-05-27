args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
require(grid)




in_f <- args[1]
out_f <- args[2]
graph <- args[3]
out_format <-args[4]
rm(args)

open_outfile <- function(out_format, out_f)
{
	if (out_format=="pdf")
	{
		pdf(out_f)
	}
	if (out_format=="png")
	{
		png(out_f)
	}
	if (out_format=="eps")
	{
		postscript(out_f)
	}
}


data=read.table(in_f, header=TRUE)
col_names=colnames(data)
n_cols=length(col_names)

labels=c()
for (i in 2:n_cols)
{
	if (col_names[i]=="aln_length")
	{
		name="alignment length"
	}
	if (col_names[i] == "identity")
	{
			name="identity (%)"
	}
	if (col_names[i]=="seq_length")
	{
			name="sequence length"
	}
	if (col_names[i]=="num_seqs")
	{
			name="# sequences"
	}
	labels=c(labels, name);
}


if (graph=="single")
{
	for (i in 2:n_cols)
	{
		full_out_f=paste(out_f,"_",col_names[i],".",out_format,sep="")
		open_outfile(out_format, full_out_f)
		plot=ggplot(data, aes(x=data[,2]))+geom_histogram(aes(fill=..count..))+scale_x_continuous(labels[i-1])
		print(plot)
		dev.off()
	}
} else {
	full_out_f=paste(out_f,".",out_format,sep="")
	open_outfile(out_format, full_out_f)
	numPlots = n_cols-1
	plotCols=1;
	if (numPlots >2)
	{
		plotCols = 2                          # Number of columns of plots
	}
	plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

	# Set up the page
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
	vplayout <- function(x, y)
		viewport(layout.pos.row = x, layout.pos.col = y)

	for (j in 2:n_cols)
	{
		full_out_f=paste(out_f,".",out_format,sep="")
		i=j-1
		# Make each plot, in the correct location
		curRow = ceiling(i/plotCols)
		curCol = (i-1) %% plotCols + 1
		this_plot=ggplot(data, aes(x=data[,j]))+geom_histogram(aes(fill=..count..))+scale_x_continuous(labels[i])
		print(this_plot, vp = vplayout(curRow, curCol ))
	}
	dev.off()
}
