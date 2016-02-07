setwd('D:\\Copy\\Bioinformatics\\2deMaster\\ComputerIntensiveBioinformatics\\cimb-project')
load("metabolite.RData")
library(reshape)
library(CMA)
library(MASS)
library(lattice)
library(ggplot2)
theme_set(theme_gray(base_size = 24))

group <- metabolites.data[,2]
metabolites <- as.matrix(metabolites.data[,-c(1:16)])
clinical <- as.matrix(metabolites.data[,c(3:16)])

# look at the expression of some metabolites
set.seed(52)
example_metabolites <- sample(1:ncol(metabolites), size = 4)
ex <- data.frame(Cancer = group, metabolites[,example_metabolites])
exdf <- melt(ex, id = c("Cancer"), variable_name = "Metabolite")
colnames(exdf)[colnames(exdf) == "value"] = "Expression"
stripplot(Expression ~ Cancer | Metabolite, exdf, grid = TRUE, group = Metabolite, 
          auto.key = TRUE, jitter.data = TRUE)

# Generate learning sets using 3-fold cross-validation with 50 repetitions
threefoldCV <- GenerateLearningsets(y=group, method="CV", fold=3, niter=50, strat=T)

# Perform feature filtering using limma ranking
metab_sel <- GeneSelection(X=metabolites, y=group, learningsets = threefoldCV, method = "limma")

show(metab_sel)
toplist(metab_sel, k = 10, iter = 1)
plot(metab_sel, iter = 1)
plot(metab_sel, iter = 2)
for (k in c(5,10,20,50)){
  sel_iter <- numeric()
  tot_iter = dim(metab_sel@rankings[[1]])[1]
  for(i in 1:tot_iter) {
    sel_iter <- c(sel_iter, toplist(metab_sel, iter = i, k = k, show = FALSE)$index) }
  print(sort(table(sel_iter), dec = T))
  barplot(sort(table(sel_iter), dec = T),las=2)
  sel_iter_df <- as.data.frame(table(sel_iter))
  print(ggplot(sel_iter_df) + geom_bar(aes(x = reorder(sel_iter, -Freq), y=Freq),stat = "identity") +
          theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ylab("Frequency") + 
          xlab('Metabolite') + ggtitle(paste("Top",k,"metabolites by LIMMA")))
}

# Analysis for top 10 metabolites
n_metab = 10

# Partial least squares: tuning
class_plsda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                              classifier = pls_ldaCMA, tuninglist = list(grids = list()))
evaluation(class_plsda)
# confusion matrix
invisible(lapply(lapply(list(class_plsda),join), ftable))
ftable(lapply(list(class_plsda),join)[[1]])
# roc curve
roc(lapply(list(class_plsda),join)[[1]])
# probability plot
plot(lapply(list(class_plsda),join)[[1]])

resultlist <- list(class_plsda, class_scda)
result <- lapply(resultlist, join)
for(i in seq(along = result)){plot(result[[i]])}

# SCDA/PAM tuning is required, but not feature selection - SLOW
tune_scda <- tune(X = metabolites, y = group, learningsets = threefoldCV, 
                  classifier = scdaCMA, grids = list())
unlist(best(tune_scda))
class_scda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                             classifier = scdaCMA, tuneres = tune_scda)
evaluation(class_scda)

# Penalized logistic regression: tuning - SLOW
class_plr <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                            classifier = plrCMA, tuninglist = list(grids = list()))
evaluation(class_plr)

# Lasso: tuning - SLOW
tune_lasso <- tune(X = metabolites, y = group, learningsets = threefoldCV, 
                   classifier = LassoCMA, grids = list())
table(unlist(best(tune_scda)))
for(i in 1:4) {plot(tune_lasso, iter = i, main = paste("iteration", i))}

class_lasso <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                              classifier = LassoCMA, tuninglist = list(grids = list()))
evaluation(class_lasso, y=group,measure = "misclassification")

# Linear, quadratic and Fisher Discriminant analyses: require feature selection, no tuning
class_fda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                            classifier = fdaCMA, genesel = metab_sel, nbgene = n_metab,comp=2)
class_lda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                            classifier = ldaCMA, genesel = metab_sel, nbgene = n_metab)
class_qda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                            classifier = qdaCMA, genesel = metab_sel, nbgene = n_metab)

# Comparing different top K's
plsda <- scda <- plr <- lasso <- fda <- lda <- qda <- list()
results <- list()
classifierlist <- list()
comparisonlist <- list()
comparisondf <- data.frame(k=numeric(), Classifier=character(), Misclassification = numeric(), Sensitivity=numeric(),
                           Specificity=numeric(), AUC=numeric(), Average_Probability=numeric())
counter <- 0
for (n_metab in seq(3,48,5)) {
  class_plsda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                                classifier = pls_ldaCMA, tuninglist = list(grids = list()))
  class_scda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                               classifier = scdaCMA, tuninglist = list(grids = list()))
  class_plr <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                              classifier = plrCMA, tuninglist = list(grids = list()))
  class_lasso <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                                classifier = LassoCMA, tuninglist = list(grids = list()))
  class_fda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                              classifier = fdaCMA, genesel = metab_sel, nbgene = n_metab,comp=2)
  class_lda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                              classifier = ldaCMA, genesel = metab_sel, nbgene = n_metab)
  class_qda <- classification(X = metabolites, y = group, learningsets = threefoldCV, 
                              classifier = qdaCMA, genesel = metab_sel, nbgene = n_metab)
  
  classifierlist[[n_metab]] <- list(class_fda,class_lda,class_qda,class_lasso,class_plr,class_plsda,class_scda)
  comparison <- compare(classifierlist[[n_metab]], plot = F, measure = c("misclassification", "sensitivity", "specificity", "auc","average probability"))
  levels(comparisondf$Classifier) <- rownames(comparison)
  for (i in 1:length(rownames(comparison))){
    comparisondf[counter+i,] <- c(as.numeric(n_metab), rownames(comparison)[i] ,c((unlist(comparison[i,]))))
  }
  counter <- counter + length(rownames(comparison))
}

# Fix dataframe: convert to numeric again
comparisondf$Misclassification <- as.numeric(comparisondf$Misclassification)
comparisondf$Sensitivity <- as.numeric(comparisondf$Sensitivity)
comparisondf$Specificity <- as.numeric(comparisondf$Specificity)
comparisondf$AUC <- as.numeric(comparisondf$AUC)
comparisondf$Average_Probability <- as.numeric(comparisondf$Average_Probability)
comparisondf$k <- as.numeric(comparisondf$k)

# Plot performance measures for different methods for different number of top k genes
ggplot(data=comparisondf) + geom_line(aes(y=Misclassification,x=k,colour=Classifier,group=Classifier),alpha=0.8)
ggplot(data=comparisondf) + geom_line(aes(y=Sensitivity,x=k,colour=Classifier,group=Classifier),alpha=0.8)
ggplot(data=comparisondf) + geom_line(aes(y=Specificity,x=k,colour=Classifier,group=Classifier),alpha=0.8)
ggplot(data=comparisondf) + geom_line(aes(y=AUC,x=k,colour=Classifier,group=Classifier),alpha=0.8)

# Compare methods for top 13 genes
comparison13 <- compare(classifierlist[[13]], plot = T, measure = c("misclassification", "sensitivity", "specificity", "auc","average probability"))
print(comparison13)
evaluation(classifierlist[[13]][[4]])
evaluation(classifierlist[[13]][[4]],measure = 'sensitivity')
evaluation(classifierlist[[13]][[4]],measure = 'specificity')
# confusion matrix
ftable(lapply(list(classifierlist[[13]][[4]]),join)[[1]])
# roc curve
roc(lapply(list(classifierlist[[13]][[4]]),join)[[1]])
# probability plot
plot(lapply(list(classifierlist[[13]][[4]]),join)[[1]])

resultlist <- list(classifierlist[[13]][[4]], class_scda)
result <- lapply(resultlist, join)
for(i in seq(along = result)){plot(result[[i]])}

# Final LASSO model
lasso_model<- LassoCMA(X = metabolites, y = group, norm.fraction = 0.1, models = T)

show(lasso_model)
ftable(lasso_model)
plot(lasso_model)

# Metabolites selected by LASSO in final model
colnames(metabolites)[lasso_test@varsel != 0]
length(colnames(metabolites)[lasso_test@varsel != 0]) # number of m
# Metabolites left out
setdiff(colnames(metabolites),colnames(metabolites)[lasso_test@varsel != 0])

### Added predictive value of metabolites compared to clinical data
library(globalboosttest)

additive_test<-globalboosttest(X=metabolites[,lasso_test@varsel !=0],Y=group,Z=clinical,nperm=25,mstop=c(100,500,1000),mstopAIC = T,plot=T,pvalueonly=FALSE)
additive_test<-globalboosttest(X=metabolites[,lasso_test@varsel !=0],Y=group,Z=model.matrix( ~ ., clin_df),nperm=25,mstop=c(100,500,1000),mstopAIC = T,plot=T,pvalueonly=FALSE)
summary(glm(as.formula(paste("Group ~ . -", paste(paste("M",1:102, sep=""), collapse= "-"))),
            data=metabolites.data[,-1],family = 'binomial'))

library(globaltest)
gt(Group~1,Group~.,data=metabolites.data[,-1],model='logistic') # clinical + metabolites vs null
gt(Group~1,as.formula(paste("Group ~ ", paste(paste("M",1:102, sep=""), collapse= "+"))),data=metabolites.data[,-1],model='logistic') # metabolites
gt(Group~1,as.formula(paste("Group ~ . -", paste(paste("M",1:102, sep=""), collapse= "-"))),data=metabolites.data[,-1],model='logistic') # clinical

a <- numeric(length(colnames(metabolites)))
a[lasso_test@varsel != 0] <- 1
which(a == 1)

gt(as.formula(paste("Group ~ ", paste(paste("M", which(a == 1) , sep=""), collapse= "+"))), # selected metabolites
   as.formula(paste("Group ~ . -", paste(paste("M", which(a == 0), sep=""), collapse= "-"))), # clinical + selected metabolites
   data=metabolites.data[,-1],model='logistic') # added value of selected metabolites

gt(as.formula(paste("Group ~ . -", paste(paste("M",1:102, sep=""), collapse= "-"))), # clinical
   as.formula(paste("Group ~ . -", paste(paste("M", which(a == 0), sep=""), collapse= "-"))), # clinical + selected metabolites
   data=metabolites.data[,-1],model='logistic') # added value of selected metabolites

# LASSO CV for clinical variables
clin_df <- metabolites.data[,3:16]
class_lasso_clinical <- classification(X =model.matrix( ~ ., clin_df), y = group, learningsets = threefoldCV, 
                              classifier = LassoCMA, tuninglist = list(grids = list()))
evaluation(class_lasso_clinical,y=group)
ftable(lapply(list(class_lasso_clinical),join)[[1]])