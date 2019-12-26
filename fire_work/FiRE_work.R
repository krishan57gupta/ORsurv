library('FiRE')
source('~/Ahuja_Lab/preprocess.R')
library(Matrix)
library(gdata)
library(reshape2)
library(ggridges)
library(ggplot2)

set.seed(0)

EXP=c("EXP0005","EXP0006","EXP0008",
      # "EXP0009","EXP0011",
      "EXP0014","EXP0017","EXP0021","EXP0023",
      "EXP0026","EXP0028",
      # "EXP0029",
      "EXP0033",
      #"EXP0034",
      "EXP0035","EXP0037","EXP0038","EXP0041",
      "EXP0050","EXP0052","EXP0053","EXP0054","EXP0056","EXP0057","EXP0064","EXP0066","EXP0071","EXP0072")
# p1=list()
# p2=list()
# p3=list()
# p4=list()
p5=list()
p6=list()
p7=list()
k=1

fire_score=list()
for(EXP_num in EXP)
{
  print(EXP_num)
true_data_1=read.csv2(paste("~/Ahuja_Lab/Files_for_Fire-20191126T105715Z-001/Files_for_Fire/",EXP_num,"_Doublets_removed.csv",sep=""),sep = ",")
true_data=true_data_1[,-1]
colnames(true_data)=colnames(true_data_1)[-1]
rownames(true_data)=true_data_1[,1]
dim(true_data)
true_data_info=read.csv2(paste("~/Ahuja_Lab/Files_for_Fire-20191126T105715Z-001/Files_for_Fire/",EXP_num,"_Final_Receptor_file.csv",sep=""),sep = ",")
true_ors_info=read.csv2(paste("~/Ahuja_Lab/Files_for_Fire-20191126T105715Z-001/Count_of_Ors/",EXP_num,"_zfpkm_or_cell_status.csv",sep=""),sep = ",")
common_cells=intersect(true_data_info[,1],colnames(true_data))
common_cells=intersect(common_cells,true_ors_info[,1])
true_cluster=true_data_info[match(common_cells,true_data_info[,1]),2]
true_ors_count=true_ors_info[match(common_cells,true_ors_info[,1]),(dim(true_ors_info)[2]-1)]
true_data=true_data[,match(common_cells,colnames(true_data))]
data=t(true_data)
data=apply(data,2,function(x) as.numeric(x))
data_mat=list("mat"=data,"gene_symbols"=colnames(data),"barcodes"=rownames(data))
preprocessedList <- ranger_preprocess(data_mat,ngenes_keep=1000)
preprocessedData <- as.matrix(preprocessedList$preprocessedData)
model <- new(FiRE::FiRE, 100, 50, 1017881, 5489, 0)
model$fit(preprocessedData)
score <- model$score(preprocessedData)
score_1=score
names(score_1)=common_cells
fire_score[[k]]=score_1

# p1[[k]]=plot(score,true_cluster,main = EXP_num)
df=data.frame("density"=score,"score"=score,"cluster"=as.factor(true_cluster),"ros_count"=true_ors_count)
df_1=data.frame("density"=score,"score"=score,"cluster"=true_cluster,"ros_count"=true_ors_count)
# p2[[k]]<-ggplot(df, aes(density, colour=cluster)) + geom_density() + ggtitle(EXP_num)
# p3[[k]]<-ggplot(df, aes(density,fill=cluster, colour=cluster)) + geom_density() + ggtitle(EXP_num)
# p4[[k]]<-ggplot(df, aes(density,fill=cluster,colour=cluster)) + geom_density(position = "stack") + ggtitle(EXP_num)
p5[[k]]<-ggplot(df, aes(density, y = cluster)) + geom_density_ridges(aes(fill = cluster)) + ggtitle(EXP_num)
p6[[k]]<-ggplot(df, aes(x=score, y = ros_count, color=cluster)) + geom_point()  + ggtitle(EXP_num)
p7[[k]]<-ggplot(df_1, aes(x=score, y = cluster, color=ros_count)) + geom_point()  + ggtitle(EXP_num)

k=k+1
}

names(fire_score)=EXP
saveRDS(fire_score,"~/Ahuja_Lab/fire_score.rds")
# pdf(paste("~/Ahuja_Lab/plot/","2D_plot.pdf",sep=""))
# p1
# dev.off()
# 
# pdf(paste("~/Ahuja_Lab/plot/","density.pdf",sep=""))
# p2
# dev.off()
# 
# pdf(paste("~/Ahuja_Lab/plot/","density_fill.pdf",sep=""))
# p3
# dev.off()
# 
# pdf(paste("~/Ahuja_Lab/plot/","density_fill_stack.pdf",sep=""))
# p4
# dev.off()

pdf(paste("~/Ahuja_Lab/plot/","density_fill_line.pdf",sep=""))
p5
dev.off()

pdf(paste("~/Ahuja_Lab/plot/","ors_counts_score_cluster.pdf",sep=""))
p6
dev.off()

pdf(paste("~/Ahuja_Lab/plot/","score_cluster_ors_count.pdf",sep=""))
p7
dev.off()




