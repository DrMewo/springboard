library(knitr)
library(reshape)
library(ggplot2)
library(gridExtra)
library(randomForest)
library(ROCR)
library(e1071)


ggplot.histogram <- function(data, title = "RNASeq Expression" )
{
	ggplot.hist <- ggplot(melt(as.matrix(data)),aes(x=value)) + geom_histogram(bins = 500) + ggtitle(title)
	return(ggplot.hist)
}

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

generate.testing.training.sets <- function(data, training.percent = .8)
{
    all.event.1 <- rownames(data$clinical)[which(data$clinical[,"event"]==1)]
    all.event.0 <- rownames(data$clinical)[which(data$clinical[,"event"]==0)]
    split.event.1 <- sample(2,length(all.event.1),replace=T,prob=c(training.percent,1-training.percent))
    split.event.0 <- sample(2,length(all.event.0),replace=T,prob=c(training.percent,1-training.percent))
    training.rows <- c(all.event.0[which(split.event.0==1)],all.event.1[which(split.event.1==1)])
    testing.rows <- c(all.event.0[which(split.event.0==2)],all.event.1[which(split.event.1==2)])
    return(list(training=list(clinical=data$clinical[training.rows,],expression=data$expression[training.rows,]),testing=list(clinical=data$clinical[testing.rows,],expression=data$expression[testing.rows,])))
}

run.classifier <- function (data, classifier = c("random forest", "svm"))
{
	data$time <- NULL
	if (classifier == "random forest")
	{
		#bestmtry <- tuneRF (data, data$event, ntreeTry=100, stepFactor=1.5, improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)
		#model <- randomForest(as.factor(event) ~ . , data = data, , mtry=bestmtry[nrow(bestmtry)-1,1], ntree=1000, keep.forest=TRUE, importance=TRUE)
		model <- randomForest(as.factor(event) ~ . , data = data, ntree=1000, keep.forest=TRUE, importance=TRUE)
	}
	if (classifier == "svm")
	{
		model <- svm(as.factor(event) ~ . , data = data)
	}
	return (model)
}

evaluate.classifier <- function (data, model)
{
	data$time <- NULL
	predict.rF <-  predict(model, type="prob", data)[,2]
	prediction.rF <- prediction (predict.rF, data$event)
	performance.rF <- performance (prediction.rF, "tpr", "fpr")
	plot(performance.rF, main= "ROC Curve for Random Forest", col=2, lwd=2)
	abline (a=0, b=1, lwd=2, lty =2, col="gray")
}

bagging.rFs <- function (data.1,data.2)
{
	data.1$time <- NULL
	data.2$time <- NULL
	control <- trainControl(method="repeatedcv", number=10, repeats=3)
	set.seed(7)
	rF.1 <- train(as.factor(event) ~ . , data = data.1, method = "rf", metric = "Accuracy", trControl = control)
	rF.2 <- train(as.factor(event) ~ . , data = data.2, method = "rf", metric = "Accuracy", trControl = control)
	bagging_results <- resamples(list(model.1=rF.1, model.2=rF.2))
	summary(bagging_results)
	dotplot(bagging_results)
}

ensemble.svm.rF <- function (predict.rF, predict.svm)
{
	predict.ensembl <- (as.numeric(as.character(predict.rF))+as.numeric(as.character(predict.svm)))/2 
}
