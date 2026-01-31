#This script integrates multiple features for gene prioritization using an extended isolation forest.

###
###
###Libraries and parameters
set.seed(1234)
options(scipen = 999)
library(data.table)
library(optparse)
library(isotree)
library(fastshap)
library(ggplot2)

###
###
###Arguments
option_list = list(
  make_option("--dataList", type="character", default=NULL,
              help="Path to file containing input data information"),

  make_option("--shapNo", type="numeric", default=NULL,
              help="Number of top predictions to compute Approximate Shapley values for"),
    
  make_option("--outName", type="character", default=NULL,
              help="Name of output")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$dataList) | is.null(opt$shapNo) | is.null(opt$outName)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###
###Functions
pfun <- function(object, newdata) {
  predict(object, newdata = newdata)
}

###
###
###Main

##
##Read and process data
inList <- fread(opt$dataList, data.table = FALSE, header = TRUE)
inList <- inList[order(inList$symEns, decreasing = TRUE),] 

mainTab <- fread(inList$path[1], header = TRUE, data.table = FALSE, select = as.numeric(unlist(strsplit(inList$columns[1],","))))
colnames(mainTab) <- unlist(strsplit(inList$colnames[1],","))

#-log10 transform
mainTab[,ncol(mainTab)] <- -log10(mainTab[,ncol(mainTab)])
for(i in 2:nrow(inList)){
	newTab <- fread(inList$path[i], header = TRUE, data.table = FALSE, select = as.numeric(unlist(strsplit(inList$columns[i],","))))
	colnames(newTab) <- unlist(strsplit(inList$colnames[i],","))
	#-log10 transform
	if(inList$trans[i] == "yes"){
		newTab[,ncol(newTab)] <- -log10(newTab[,ncol(newTab)])
	}
	mainTab <- merge(mainTab, newTab, by = c(unlist(strsplit(inList$symEns[i],","))), all.x = TRUE, all.y = TRUE)
}

#Impute missing to median
mainTab <- mainTab[!(is.na(mainTab$ensemblID)),]
for(i in 3:ncol(mainTab)){
	mainTab[,i][is.na(mainTab[,i])] <- median(mainTab[,i], na.rm = TRUE)
}

##
##IF
trainMod <- isolation.forest(mainTab[,3:ncol(mainTab)], ndim=2, ntrees=100, nthreads = 5)
save(trainMod, file = paste0("../data/out/",opt$outName,".RData"))
mainTab$EIF <- predict(trainMod, mainTab)
mainTab$rank <- rank(-mainTab$EIF)
write.table(mainTab, paste0("../data/out/",opt$outName,"InclFeatures.txt"), sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = mainTab$ensemblID, EIF = mainTab$EIF), paste0("../data/out/",opt$outName,".txt"), sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE, quote = FALSE)

##
##Shap values for top
mainTab <- mainTab[order(-mainTab$EIF),]
ex <- data.frame(symbol = c(), Feature = c(), shap = c())
for(i in 1:opt$shapNo){
  (exShort <- fastshap::explain(trainMod, X = mainTab[,c(3:(ncol(mainTab) -2))] , newdata = mainTab[i,c(3:(ncol(mainTab) -2))], pred_wrapper = pfun, 
                           adjust = TRUE, nsim = 1000))
  symbol <- rep(mainTab$symbol[i], ncol(data.frame(exShort)))
  featShort <- names(data.frame(exShort))
  shapShort <- as.numeric(exShort)
  exShort <- data.frame(symbol = symbol, Feature = featShort,shap = shapShort)
  ex <- rbind(ex,exShort)
}

##
##Write output
write.table(ex, paste0("../data/out/",opt$outName,"Shap",i,".txt"), sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE, quote = FALSE)

##
##Plot feature contributions
levelOrder <- rev(unique(ex$symbol))
p <- ggplot(ex, aes(fill = Feature, x = shap, y = factor(symbol, level = levelOrder))) +
  geom_bar(position="stack", stat="identity")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("") +
  xlab("Approximate Shapley value") +
  theme_classic()+
  theme(axis.text.y = element_text(size = rel(0.8)))+
  theme(axis.text.y = element_text(face = "italic"))
ggsave(file = paste0("../data/out/",opt$outName,"Shap",i,".pdf"))
ggsave(file = paste0("../data/out/",opt$outName,"Shap",i,".png"), width = 8, height = 5, dpi = 600)
