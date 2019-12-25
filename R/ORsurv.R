##' @title survival of patients using ors on TCGA data
##' @description survival of patients using ors on TCGA data
##' @param folder The normalised and filtered
##' @param matrix_1 The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param matrix_2 Labels for the two sub-populations
##' @param matrix_3 optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param t1 Labels for the two sub-populations
##' @param t2 optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param t3 The number of cores to be used
##' @param name1 The number of cores to be used
##' @param name2 Labels for the two sub-populations
##' @param name3 optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param survival_info The number of cores to be used
##' @param method The number of cores to be used
##' @param month_limit Labels for the two sub-populations
##' @param p_limit optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param selected_label The number of cores to be used
##' @return surv_plot pValues A vector containing FDR adjusted p significance values
##' @export surv_med
ORsurv<-function(folder,matrix_1,matrix_2,matrix_3,t1=-3,t2=3,t3=3,name1="EXP",name2="EXP_TCGA",name3="Survival_plot",survival_info,method,month_limit,p_limit=2,selected_label)
{
  data=matrix_filter(matrix_1=matrix_1,matrix_2=matrix_2,matrix_3=matrix_3,survival_info=survival_info)
  EXP=num2bin(mat=data$EXP_ZFPKM,t=t1)
  EXP_clusters=cluster_label(folder=folder,mat =EXP,t=t2,name=name1)
  EXP_sig_mat=signature_matrix(mat=data$EXP_TPM,label=EXP_clusters)
  corr_mat=distance_matrix(mat_1=EXP_sig_mat,mat_2=data$TCGA_TPM,method=method)
  surv_med=survival_plot(folder=folder,mat=corr_mat,mat_2=data$TCGA_TPM,s=data$survival_info,t=t3,name2=name2,name3=name3,month_limit=month_limit,p_limit=p_limit,selected_label=selected_label)
  return(surv_med)
}
##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param B optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param A optional, 10 as more accuracy
##' @return cos_sim pValues A vector containing FDR adjusted p significance values
cosine_formula<-function(A,B)
{
  cos_sim=sum(A*B)/sqrt(sum(A^2)*sum(B^2))
  return(cos_sim)
}
##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param labels optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @return new_color pValues A vector containing FDR adjusted p significance values
color_map<-function(labels)
{
  all_color=list("1"="red","2"="green","3"="blue","4"="black","5"="yellow","6"="orange","7"="purple","8"="violet","9"="gold","10"="pink")
  new_color=all_color[seq_len(length(labels))]
  names(new_color)=labels
  new_color=list("map"=new_color,"color"=unlist(new_color),"labels"=names(new_color))
  return(new_color)
}
##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param matrix_1 The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param matrix_2 Labels for the two sub-populations
##' @param matrix_3 Labels for the two sub-populations
##' @param survival_info optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @return filtered_matrix pValues A vector containing FDR adjusted p significance values
matrix_filter<-function(matrix_1,matrix_2,matrix_3,survival_info)
{
  common_ors=intersect(rownames(matrix_1),rownames(matrix_2))
  common_ors=intersect(common_ors,rownames(matrix_3))
  common_ors=sort(common_ors)
  common_sample_2=intersect(colnames(matrix_3),rownames(survival_info))
  survival_info=survival_info[match(common_sample_2,rownames(survival_info)),]
  surv_r=rownames(survival_info)
  surv_c=colnames(survival_info)
  survival_info=as.matrix(survival_info)
  survival_info=apply(survival_info,2,as.numeric)
  rownames(survival_info)=surv_r
  colnames(survival_info)=surv_c
  r=rownames(matrix_1)
  c=colnames(matrix_1)
  matrix_1=apply(matrix_1,2,as.numeric)
  rownames(matrix_1)=r
  colnames(matrix_1)=c
  r=rownames(matrix_2)
  c=colnames(matrix_2)
  matrix_2=apply(matrix_2,2,as.numeric)
  rownames(matrix_2)=r
  colnames(matrix_2)=c
  r=rownames(matrix_3)
  c=colnames(matrix_3)
  matrix_3=apply(matrix_3,2,as.numeric)
  rownames(matrix_3)=r
  colnames(matrix_3)=c
  matrix_3=matrix_3[match(common_ors,rownames(matrix_3)),match(common_sample_2,colnames(matrix_3))]
  matrix_3=matrix_3[,apply(matrix_3,2,function(x) {
    if(stats::sd(x)==0)
      return(FALSE)
    else
      return(TRUE)
  })]
  common_sample_2=intersect(colnames(matrix_3),rownames(survival_info))
  matrix_3=matrix_3[match(common_ors,rownames(matrix_3)),match(common_sample_2,colnames(matrix_3))]
  survival_info=survival_info[match(common_sample_2,rownames(survival_info)),]
  filtered_matrix=list("EXP_ZFPKM"=matrix_1[match(common_ors,rownames(matrix_1)),],
                       "EXP_TPM"=log2(matrix_2[match(common_ors,rownames(matrix_2)),]+1),
                       "TCGA_TPM"=log2(matrix_3+1),
                       "survival_info"=survival_info)
  return(filtered_matrix)
}

##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param mat The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param t Labels for the two sub-populations
##' @return label pValues A vector containing FDR adjusted p significance values
num2bin<-function(mat,t=-3) # here t should be less then 1
{
  mat_1=mat
  mat=apply(mat_1,2,as.numeric)
  colnames(mat)<-colnames(mat_1)
  rownames(mat)<-rownames(mat_1)
  mat[mat>t]=1
  mat[mat<=t]=0
  return(mat)
}

##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param folder The normalised and filtered, read count matrix, with row names
##' @param mat The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param t Labels for the two sub-populations
##' @param name Labels for the two sub-populations
##' @param label_2 optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @return label pValues A vector containing FDR adjusted p significance values
cluster_label<-function(folder,mat,t=3,name="EXP",label_2)
{

  heatmap_info=pheatmap::pheatmap(mat,cluster_rows=F)
  label=stats::cutree(heatmap_info$tree_col, k=t)


  sig_mat=signature_matrix(mat = mat,label = label)
  heatmap_info=NMF::aheatmap(sig_mat,file=paste(folder,name,"_binary_probability.pdf",sep=""),Colv = F,Rowv=F,txt = round(sig_mat,2))

  # p=prcomp(mat)
  # scatter_data = data.frame("PC1"=p$rotation[,1],"PC2"=p$rotation[,2],"labels"=factor(label))
  # g=ggplot2::ggplot(scatter_data, ggplot2::aes(x=PC1, y=PC2,color=labels)) + geom_point(size=.5)
  # grDevices::pdf(paste(folder,name,"_binary_scatter",".pdf",sep=""))
  # g
  # grDevices::dev.off()

  col_annotation = data.frame(label_1 = factor(label))
  rownames(col_annotation)<-colnames(mat)
  ann_colors=list("label_1"=color_map(unique(label))$color)
  heatmap_info=NMF::aheatmap(mat,annCol = col_annotation,labCol = "",color = factor(c("grey","#556b2f")),file=paste(folder,name,"_binary.pdf",sep=""), annColors=ann_colors)
  return(label)
}
##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param mat The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param label for the two sub-populations
##' @return sigmat pValues A vector containing FDR adjusted p significance values
signature_matrix<-function(mat,label)
{
  unique_label<-unique(label)
  sig_mat=matrix(0,ncol=length(unique_label),nrow=dim(mat)[1])
  sig_mat_name=c()
  for(i in seq_len(length(unique_label)))
  {
    sig_mat_name=c(sig_mat_name,paste("C",unique_label[i],sep="_"))
    sub_mat=mat[,which(label%in%unique_label[i])]
    if (is.vector(sub_mat))
      sig_mat[,i]=sub_mat
    else
      sig_mat[,i]=apply(sub_mat,1,mean)
  }
  colnames(sig_mat)=sig_mat_name
  rownames(sig_mat)=rownames(mat)
  return(sig_mat)
}

##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param mat_1 The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param mat_2 Labels for the two sub-populations
##' @param method optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @return mat_3 pValues A vector containing FDR adjusted p significance values
distance_matrix<-function(mat_1,mat_2, method="cosine")
{
  mat_3=matrix(0,nrow = dim(mat_1)[2],ncol=dim(mat_2)[2])
  if(method=="cosine")
  {
    for(i in seq_len(dim(mat_1)[2]))
    {
      for(j in seq_len(dim(mat_2)[2]))
      {
        mat_3[i,j]=cosine_formula(mat_1[,i],mat_2[,j])
      }
    }
  }
  if(method=="hamming")
  {
    for(i in seq_len(dim(mat_1)[2]))
    {
      for(j in seq_len(dim(mat_2)[2]))
      {
        mat_3[i,j]=e1071::hamming.distance(mat_1[,i],mat_2[,j])
      }
    }
  }
  if(method=="cor")
  {
    for(i in seq_len(dim(mat_1)[2]))
    {
      for(j in seq_len(dim(mat_2)[2]))
      {
        mat_3[i,j]=stats::cor(mat_1[,i],mat_2[,j])
      }
    }
  }
  rownames(mat_3)=colnames(mat_1)
  colnames(mat_3)=colnames(mat_2)
  return(mat_3)
}

##' @title A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param folder The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param mat The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param mat_2 Labels for the two sub-populations
##' @param s optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param t Labels for the two sub-populations
##' @param name2 Labels for the two sub-populations
##' @param name3 optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param month_limit Labels for the two sub-populations
##' @param p_limit optional, 10 as more accuracy with little slow (calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param selected_label The number of cores to be used
##' @return surv_plot pValues A vector containing FDR adjusted p significance values
##' @export surv_med
survival_plot<-function(folder,mat,mat_2,s,t=3,name2="EXP_TCGA",name3="Survival_plot",month_limit,p_limit, selected_label)
{
  ss=s
  ss[,1]=ss[,1]/30
  heatmap_info=pheatmap::pheatmap(mat,cluster_rows=F)
  label=stats::cutree(heatmap_info$tree_col, k=t)

  # col_annotation = data.frame(label = factor(label))
  # rownames(col_annotation)<-colnames(mat)
  ann_colors=list("label"=color_map(unique(label))$color)
  # heatmap_info=NMF::aheatmap(mat,annCol = col_annotation,labCol = "",file=paste(folder,name2,".pdf",sep=""),annColors = ann_colors)

  survival_mat=data.frame(t(mat),time=ss[,1],status=ss[,2],label=label)

  # avg_data=signature_matrix(mat,label)
  # heatmap_info=NMF::aheatmap(avg_data,Rowv = F,Colv = F,file=paste(folder,name2,"_avg_on_projection",".pdf",sep=""),txt=round(avg_data,2))
  #
  # new_l=c()
  # for(i in seq_len(dim(avg_data)[2]))
  #   new_l=c(new_l,rep(colnames(avg_data)[i],dim(avg_data)[1]))
  # avg_data_select_df=data.frame("x"=c(avg_data),label=factor(new_l))
  # p<-ggplot2::ggplot(avg_data_select_df, ggplot2::aes(x=label, y="x",fill=label)) +
  #   ggplot2::geom_boxplot(position=ggplot2::position_dodge(1))+ggplot2::geom_dotplot(binaxis='y', stackdir='center', binwidth=.001)
  # grDevices::pdf(paste(folder,name3,"_boxplot_on_projection",".pdf",sep=""))
  # p
  # grDevices::dev.off()
  #
  # p=prcomp(mat)
  # avg_data_select_df2=data.frame("PC1"=p$rotation[,1],"PC2"=p$rotation[,2],"label"=factor(label))
  # g=ggplot2::ggplot(avg_data_select_df2, ggplot2::aes(x=PC1, y=PC2,color=label)) + geom_point(size=.5)
  # grDevices::pdf(paste(folder,name3,"_scatter_on_projection",".pdf",sep=""))
  # g
  # grDevices::dev.off()
  #
  # avg_data=signature_matrix(mat_2,label)
  # heatmap_info=NMF::aheatmap(avg_data,Rowv = F,Colv = F,file=paste(folder,name2,"_avg_on_TCGA",".pdf",sep=""),txt=round(avg_data,2))
  #
  # new_l=c()
  # for(i in seq_len(dim(avg_data)[2]))
  #   new_l=c(new_l,rep(colnames(avg_data)[i],dim(avg_data)[1]))
  # avg_data_select_df=data.frame("x"=c(avg_data),label=factor(new_l))
  # p<-ggplot2::ggplot(avg_data_select_df, ggplot2::aes(x=label, y="x",fill=label)) +
  #   ggplot2::geom_boxplot(position=ggplot2::position_dodge(1))+ggplot2::geom_dotplot(binaxis='y', stackdir='center', binwidth=.001)
  # grDevices::pdf(paste(folder,name3,"_boxplot_on_TCGA",".pdf",sep=""))
  # p
  # grDevices::dev.off()
  #
  # p=prcomp(mat_2)
  # avg_data_select_df2=data.frame("PC1"=p$rotation[,1],"PC2"=p$rotation[,2],"label"=factor(label))
  # g=ggplot2::ggplot(avg_data_select_df2, ggplot2::aes(x=PC1, y=PC2,color=label)) + geom_point(size=.5)
  # grDevices::pdf(paste(folder,name3,"_scatter_on_TCGA",".pdf",sep=""))
  # g
  # grDevices::dev.off()

  fit<-survival::survfit(survival::Surv(time,status) ~ label,data=survival_mat)
  pair_pval=survminer::pairwise_survdiff(survival::Surv(time,status) ~ label,data=survival_mat)
  message=""
  message_last="Strata   "
  dark_labels=c()
  for(jj in seq_len(dim(pair_pval$p.value)[1]))
  {
    message=paste(message,"\n","    ",rownames(pair_pval$p.value)[jj],"       ",sep="")
    for(kk in seq_len(dim(pair_pval$p.value)[2]))
    {
      if(jj>=kk)
      {
        if(pair_pval$p.value[jj,kk]<=p_limit)
        {
          dark_labels=c(dark_labels,colnames(pair_pval$p.value)[kk])
          dark_labels=c(dark_labels,rownames(pair_pval$p.value)[jj])
        }
        if(stringr::str_length(paste(round(pair_pval$p.value[jj,kk],2)))==4)
          message=paste(message," ",round(pair_pval$p.value[jj,kk],2)," ",sep="")
        if(stringr::str_length(paste(round(pair_pval$p.value[jj,kk],2)))==3)
          message=paste(message," ",round(pair_pval$p.value[jj,kk],2),"   ",sep="")
        if(stringr::str_length(paste(round(pair_pval$p.value[jj,kk],2)))==1)
          message=paste(message," ",round(pair_pval$p.value[jj,kk],2),"      ",sep="")
      }
    }
    message_last=paste(message_last," ",colnames(pair_pval$p.value)[jj],"       ",sep="")
  }
  message=paste(message,"\n",message_last,sep="")
  g=survminer::ggsurvplot(fit, data = survival_mat ,pval = TRUE, risk.table = TRUE, pval.coord=c(0,.3),censor=FALSE
               ,palette = ann_colors$label, legend.labs =names(ann_colors$label),surv.median.line="hv",size=.5)
  surv_med=survminer::surv_median(fit, combine = FALSE)
  grDevices::pdf(paste(folder,name3,".pdf",sep=""))
  g
  grDevices::dev.off()


  dark_labels_unique=unique(dark_labels) # active when thersold on the basis of p values
  if(p_limit>1)
  dark_labels_unique=selected_label # active when thersold on the basis of selected values

  mat_new=mat[,which(label%in%dark_labels_unique)]
  mat_2_new=mat_2[,which(label%in%dark_labels_unique)]
  ss_new=ss[which(label%in%dark_labels_unique),]
  label_new=label[which(label%in%dark_labels_unique)]

  mat=mat_new[,ss_new[,1]<month_limit]
  mat_2=mat_2_new[,ss_new[,1]<month_limit]
  ss=ss_new[  ss_new[,1]<month_limit,]
  label=label_new[  ss_new[,1]<month_limit]


  col_annotation = data.frame(label = factor(label))
  rownames(col_annotation)<-colnames(mat_2)
  ann_colors=list("label"=color_map(unique(label))$color)
  heatmap_info=NMF::aheatmap(mat,annCol = col_annotation,labCol = "",file=paste(folder,"filtered",name2,".pdf",sep=""),annColors = ann_colors)
  survival_mat=data.frame(t(mat),time=ss[,1],status=ss[,2],label=label)

  avg_data=signature_matrix(mat,label)
  heatmap_info=NMF::aheatmap(avg_data,Rowv = F,Colv = F,file=paste(folder,"filtered",name2,"_avg_on_projection",".pdf",sep=""),txt=round(avg_data,2))

  new_l=c()
  for(i in seq_len(dim(avg_data)[2]))
    new_l=c(new_l,rep(colnames(avg_data)[i],dim(avg_data)[1]))
  avg_data_select_df=data.frame("x"=c(avg_data),label=factor(new_l))
  p<-ggplot2::ggplot(avg_data_select_df, ggplot2::aes(x=label, y="x",fill=label)) +
    ggplot2::geom_boxplot(position=ggplot2::position_dodge(1))+ggplot2::geom_dotplot(binaxis='y', stackdir='center', binwidth=.001)
  grDevices::pdf(paste(folder,"filtered",name3,"_boxplot_on_projection",".pdf",sep=""))
  p
  grDevices::dev.off()

  # p=prcomp(mat)
  # avg_data_select_df2=data.frame("PC1"=p$rotation[,1],"PC2"=p$rotation[,2],"label"=factor(label))
  # g=ggplot2::ggplot(avg_data_select_df2, ggplot2::aes(x=PC1, y=PC2,color=label)) + geom_point(size=.5)
  # grDevices::pdf(paste(folder,"filtered",name3,"_scatter_on_projection",".pdf",sep=""))
  # g
  # grDevices::dev.off()

  avg_data=signature_matrix(mat_2,label)
  heatmap_info=NMF::aheatmap(avg_data,Rowv = F,Colv = F,file=paste(folder,"filtered",name2,"_avg_on_TCGA",".pdf",sep=""),txt=round(avg_data,2))

  new_l=c()
  for(i in seq_len(dim(avg_data)[2]))
    new_l=c(new_l,rep(colnames(avg_data)[i],dim(avg_data)[1]))
  avg_data_select_df=data.frame("x"=c(avg_data),label=factor(new_l))
  p<-ggplot2::ggplot(avg_data_select_df, ggplot2::aes(x=label, y="x",fill=label)) +
    ggplot2::geom_boxplot(position=ggplot2::position_dodge(1))+ggplot2::geom_dotplot(binaxis='y', stackdir='center', binwidth=.001)
  grDevices::pdf(paste(folder,"filtered",name3,"_boxplot_on_TCGA",".pdf",sep=""))
  p
  grDevices::dev.off()

  # p=prcomp(mat_2)
  # avg_data_select_df2=data.frame("PC1"=p$rotation[,1],"PC2"=p$rotation[,2],"label"=factor(label))
  # g=ggplot2::ggplot(avg_data_select_df2, ggplot2::aes(x=PC1, y=PC2,color=label)) + geom_point(size=.5)
  # grDevices::pdf(paste(folder,"filtered",name3,"_scatter_on_TCGA",".pdf",sep=""))
  # g
  # grDevices::dev.off()

  survival_mat=data.frame(t(mat),time=ss[,1],status=ss[,2],label=label)
  fit<-survival::survfit(survival::Surv(time,status) ~ label,data=survival_mat)
  pair_pval=survminer::pairwise_survdiff(survival::Surv(time,status) ~ label,data=survival_mat)
  message=""
  message_last="Strata   "
  for(jj in seq_len(dim(pair_pval$p.value)[1]))
  {
    message=paste(message,"\n","    ",rownames(pair_pval$p.value)[jj],"       ",sep="")
    for(kk in seq_len(dim(pair_pval$p.value)[2]))
    {
      if(jj>=kk)
      {
        if(stringr::str_length(paste(round(pair_pval$p.value[jj,kk],2)))==4)
          message=paste(message," ",round(pair_pval$p.value[jj,kk],2)," ",sep="")
        if(stringr::str_length(paste(round(pair_pval$p.value[jj,kk],2)))==3)
          message=paste(message," ",round(pair_pval$p.value[jj,kk],2),"   ",sep="")
        if(stringr::str_length(paste(round(pair_pval$p.value[jj,kk],2)))==1)
          message=paste(message," ",round(pair_pval$p.value[jj,kk],2),"      ",sep="")
      }
    }
    message_last=paste(message_last," ",colnames(pair_pval$p.value)[jj],"       ",sep="")
  }
  message=paste(message,"\n",message_last,sep="")
  g=survminer::ggsurvplot(fit, data = survival_mat ,pval = TRUE, risk.table = TRUE, pval.coord=c(0,.3),censor=FALSE
               ,palette = color_map(unique(label))$color, legend.labs = names(ann_colors$label),surv.median.line="hv",size=.5)
  grDevices::pdf(paste(folder,"filtered",name3,".pdf",sep=""))
  g
  grDevices::dev.off()
  return(surv_med)
}
