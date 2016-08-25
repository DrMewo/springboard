library(knitr)
library(reshape)



normalize.expression.data <- function(data)
{
    data[is.na(data)] <- 0
    tmp <- data-apply(data,2,median)
    normalized.data <- tmp - apply(tmp,1,median)
    return(normalized.data)
}



select.variable.genes <- function (data, quantile.min = .33, z.min = .5)
{
	sd.by.gene <- apply(data,2,sd)
	sd.min <- as.numeric(quantile(sd.by.gene,quantile.min))
	mean.by.gene <- apply(data,2,mean)
	z.score <- (mean.by.gene-mean(mean.by.gene))/sd(mean.by.gene)
	genes <- colnames(data)[which(abs(z.score) >= z.min & sd.by.gene >= sd.min)]
	genes.rank <- rank(abs(z.score[genes]))
	return (list(genes = genes, genes.rank = genes.rank))
}



remove.correlated.genes.1 <- function (data, correlation.max = .5)
{
	cor.matrix <- cor(data$expression.data)
	cor.matrix[upper.tri(cor.matrix,diag=TRUE)]<-NA
	cor.table <- setNames(melt(cor.matrix), c('gene_1', 'gene_2', 'value'))
	cor.table[which(abs(cor.table$value) <= correlation.max),] <- NA
	cor.table<-cor.table[which(!is.na(cor.table$value)),]
	if (dim(cor.table)[1] > 0)
	{
		genes <- remove.correlated.genes.2(list(cor.table = cor.table, genes.rank = data$genes.rank, all.genes = colnames(data$expression.data)))
	}
	else
	{
		genes <- colnames(cor.matrix)
	}
	return (genes)
}

remove.correlated.genes.2 <- function (data.2)
{
	cor.table <- data.2$cor.table
	cor.table <- cor.table[order(-abs(cor.table$value)),]
	cor.table$rank_1 <- data.2$genes.rank[cor.table$gene_1]
	cor.table$rank_2 <- data.2$genes.rank[cor.table$gene_2]
	cor.table$remove <- rep('NA',dim(cor.table)[1])
	cor.table$remove <- ifelse((cor.table$rank_1 >= cor.table$rank_2),as.character(cor.table$gene_1),as.character(cor.table$gene_2))
	genes<-data.2$all.genes[which(!(data.2$all.genes%in%unique(cor.table$remove)))]
	return(genes)
}
