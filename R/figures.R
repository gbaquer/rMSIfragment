#figures.R : Functions to create figures
#ROC curve
fig_dir<-"/home/gbaquer/msidata/1. In-source Fragmentation/1.6. Paper/Figures/"
td_roc<-function(scores,weights,mz){
  r<-list()
  r$scores<-sort(unique(scores))

  r$curve<-data.frame(fpr=rep(0,length(r$scores)))
  r$curve$fpr<-sapply(r$scores,function(x)sum(!weights[scores>=x])/sum(!weights))
  r$curve$tpr<-sapply(r$scores,function(x)sum(weights[scores>=x])/sum(weights))

  r$curve2<-data.frame(hits=rep(0,length(r$scores)))
  r$curve2$hits<-sapply(r$scores,function(x)sum(scores>=x)/length(unique(mz)))
  r$curve2$fdr<-sapply(r$scores,function(x)sum(!weights[scores>=x])/sum(weights[scores>=x]))
  r$curve2<-r$curve2[r$curve2$hits>0.05,]
  i<-order(r$curve$fpr)
  r$auc.integral<-sum(diff(r$curve$fpr[i])*(head(r$curve$tpr[i],-1)+tail(r$curve$tpr[i],-1)))/2
  return(r)
}
plot_roc <- function(res,metrics,names=metrics,mode="ROC"){
  if( mode=="ROC")
    roc<-lapply(metrics,function(x)with(res,td_roc(eval(parse(text=x)),db=="T",mz_i)))
  else
    roc<-lapply(metrics,function(x)with(res,PRROC::pr.curve(scores.class0 = eval(parse(text=x)),weights.class0 = db=="T",curve=T)))

  if(mode=="ROC"){
    plot(0,type="l",col=1,xlab="Hits per mz",ylab="False Discovery Rate (Decoy/Target)",xlim=c(0,18),ylim=c(0,1))
    for (i in 1:length(roc)){
      lines(roc[[i]]$curve2,col=i,lw=2)
    }
    legend(8.5,1.05,paste(names,round(sapply(roc,function(x)x$auc.integral),2)),col=1:(length(roc)),lw=2)
  }

  plot(0:1,0:1,type="l",col=1,asp=1,ylab="True Positive Rate (Target Rate)",xlab="False Positive Rate (Decoy Rate)")
  for (i in 1:length(roc)){
    lines(roc[[i]]$curve,col=i,lw=2)
  }
  legend(0.5,0.4,paste(names,round(sapply(roc,function(x)x$auc.integral),2)),col=1:(length(roc)),lw=2)


}
plot_figure1<-function(res,metrics,names=metrics,label=""){
  #Target Decoy Validation
  roc<-lapply(metrics,function(x)with(res,td_roc(eval(parse(text=x)),db=="T",mz_i)))
  pr<-lapply(metrics,function(x)with(res,PRROC::pr.curve(scores.class0 = eval(parse(text=x)),weights.class0 = db=="T",curve=T)))

  df<-data.frame(x=unlist(sapply(roc,function(x)x$curve$fpr)),y=unlist(sapply(roc,function(x)x$curve$tpr)),var=rep(names,times=sapply(roc,function(x)nrow(x$curve))))

  pa<-ggplot(df,aes(x,y,col=var))+geom_line(size=2)+
    xlab("False Positive Rate (Decoy Rate)")+ylab("True Positive Rate (Target Rate)")+
    theme_bw()+theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("A")+
    scale_color_discrete(name = "Score", breaks = names, labels = paste(names," (",round(sapply(roc,function(x)x$auc.integral),2)," AUC)",sep = ""))

  df<-data.frame(x=unlist(sapply(pr,function(x)x$curve[,1])),y=unlist(sapply(pr,function(x)x$curve[,2])),var=rep(names,times=sapply(pr,function(x)nrow(x$curve))))

  pb<-ggplot(df,aes(x,y,col=var))+geom_line(size=2)+
    xlab("Recall")+ylab("Precision")+
    theme_bw()+theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("B")+
    scale_color_discrete(name = "Score", breaks = names,labels = paste(names," (",round(sapply(pr,function(x)x$auc.integral),2)," AUC)",sep = ""))



  df<-data.frame(x=unlist(sapply(roc,function(x)x$curve2$hits)),y=unlist(sapply(roc,function(x)x$curve2$fdr)),var=rep(names,times=sapply(roc,function(x)nrow(x$curve2))))

  pc<-ggplot(df,aes(x,y,col=var))+geom_line(size=2)+
    xlab("Annotations per MS signal")+ylab("False Discovery Rate (Decoy/Target)")+
    theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("C")+
    scale_color_discrete(name = "Score", breaks = names,labels = paste(names," (",round(sapply(roc,function(x)x$auc.integral),2)," AUC)",sep = ""))

  p=grid.arrange(pa,pc,nrow=1)
  ggsave(paste("Fig1_Target_Decoy_Validation",".tiff",sep=label),p,path=fig_dir,width=20,height=10)
}

plot_figure2<-function(r,t,s="S",steps=1000,label=""){
  #Manual validation
  df<-data.frame(o=rep(0,steps),o_hplc=rep(0,steps),th=seq(min(r[,s]),max(r[,s]),length.out=steps))
  df$n<-sapply(df$th,function(x)sum(r[,s]>=x))
  for(i in seq_along(df$th)){
    m<-manual_validation(r[r[,s]>=df$th[i],],t)
    df$o[i]<-sum(m$table$match)/nrow(m$table)
    df$o_hplc[i]<-sum(subset(m$table,HPLC_adduct!="_")$match)/nrow(subset(m$table,HPLC_adduct!="_"))
  }
  c<-max(df$n)
  dfp<-data.frame(threshold=rep(df$th,3),y=with(df,c(o,o_hplc,n/c)),var=rep(c("Matches (%)","Matches(%) (HPLC validated)","Annotations per MS signal"),each=steps))
  pa<-ggplot(dfp,aes(threshold,y,col=var))+geom_line(size=2)+
    scale_y_continuous(
      name = "Matches to Garate et al. 2020 (%)",
      sec.axis = sec_axis(~.*c/length(pks$mass), name="Annotations per MS signal")
    ) +xlab("Ranking Score (S) threshold")+theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("A")

  m<-manual_validation(r,t)
  dfp<-data.frame(score=m$res[,s],type=as.factor(match(m$res$match_HPLC_adduct,c("","_"),nomatch=3)))
  pb<-ggplot(dfp,aes(type,score,col=type))+geom_violin()+geom_point()+
    stat_summary(fun = "median",
                 geom = "crossbar",
                 width = 0.5)+
    xlab("")+ylab("Ranking Score (S)")+theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("B")

  p=grid.arrange(pa,pb)
  ggsave(paste("Fig2_Manual_Validation",".tiff",sep=label),p,path=fig_dir,width=10,height=10)
}

require(data.table)
plot_figure2_topN<-function(r,t,s="S",steps=5,label=""){
  #Manual validation
  df<-data.frame(o=rep(0,steps),o_hplc=rep(0,steps),th=steps:1)
  r<-manual_validation(r,t)$res
  for(i in seq_along(df$th)){
    m<-manual_validation(topN(r,s,df$th[i]),t)
    df$o[i]<-sum(m$table$match)/nrow(m$table)
    df$o_hplc[i]<-sum(subset(m$table,HPLC_adduct!="_")$match)/nrow(subset(m$table,HPLC_adduct!="_"))
  }
  dfp<-data.frame(threshold=rep(df$th,2),y=with(df,c(o,o_hplc)),var=rep(c("Matches (%)","Matches(%) (HPLC validated)"),each=steps))
  pa<-ggplot(dfp,aes(threshold,y,col=var))+geom_line(size=2)+geom_point(size=10)+xlab("Top N for each mz")+ylab("Matches to Garate et al. 2020 (%)")+theme_bw()+theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                                                                                                                                      panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("A")+ylim(0,0.75)

  m<-manual_validation(r,t)
  dfp<-data.frame(score=m$res[,s],type=as.factor(match(m$res$match_HPLC_adduct,c("","_"),nomatch=3)))
  pb<-ggplot(dfp,aes(type,score,col=type))+geom_violin()+geom_point()+
    stat_summary(fun = "median",
                 geom = "crossbar",
                 width = 0.5)+
    xlab("")+ylab("Ranking Score (S)")+theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("B")

  p=pa
  ggsave(paste("Fig2_Manual_Validation_TopN",".tiff",sep=label),p,path=fig_dir,width=10,height=10)
}

plot_figure3<-function(pks,res,mz,topN=NA,score="lipid_occurences*(1+correlation)",label=""){
  res<-res[order(with(res,eval(parse(text=score))),decreasing=T),]
  i<-which.min(abs(mz-pks$mass))
  x<-subset(res,mz_i==i)
  if(!is.na(topN)&nrow(x)>0){
    x<-x[1:min(nrow(x),topN),]
  }
  k<-list()
  for(i in 1:nrow(x)){
    y<-subset(res,lipid_id==x$lipid_id[i])
    z<-matrix(0,nrow(y)+1,nrow(y)+1)
    n=with(rbind(x[i,],y),paste(round(experimental_mz,2),"\n[",adduct," ",fragmentation,"]",sep=""))
    n[1]<- paste(x$abbreviation[i],x$formula[i],n[1],paste("LO:",x$lipid_occurences[i],"C:",round(x$correlation[i],2)),sep="\n")
    z[1,-1]=1
    k[[i]]<-list(z=z,n=n,t=c(0,1+as.numeric(y$fragmentation=="")))
  }
  net = network(do.call(magic::adiag,sapply(k,function(x)x$z)), directed = T)
  pa<-ggnet2(net,mode="kamadakawai",label=unlist(sapply(k,function(x)x$n)),shape=15,label.size=3.5,size=20,color=1+unlist(sapply(k,function(x)x$t)))+
    ggtitle("A")+
    scale_color_brewer(palette="Set3",labels = c("Annotation for current mz","Fragment","Parental ion"))



  y<-rbind(x[1,],subset(res,lipid_id==x$lipid_id[1]))
  pb<-grid.arrange(top=grid::textGrob("B", x = 0, hjust = 0),grobs=lapply(1:nrow(y),function(i)ggplot_peak_image(one_pks(pks,4),y$mz_i[i],k[[1]]$n[i],k[[1]]$t[i]+1)))

  p<-grid.arrange(pa,pb,ncol=1)
  ggsave(paste("Fig3_Example_Annotations",".tiff",sep=label),p,path=fig_dir,width=10,height=20)
}
ggplot_mean_spectra<- function(pks,i1=NA,i2=NA,i3=NA){
  col=rep("black",length(pks$mass))

  p=ggplot(data.frame(mz=pks$mass,i=apply(pks$intensity,2,mean)),
           aes(x=mz, ymax=i, ymin=0,col=col)) +
    geom_linerange()

  if(!is.na(i1))
    p<-p+ geom_vline(xintercept = pks$mass[as.numeric(i1)], linetype="dotted",
                     color = "red", size=1.5)
  if(!is.na(i2))
    p<-p+ geom_vline(xintercept = pks$mass[as.numeric(i2)], linetype="dotted",
                     color = "green", size=1.5)
  if(!is.na(i3))
    p<-p+ geom_vline(xintercept = pks$mass[as.numeric(i3)], linetype="dotted",
                     color = "blue", size=1.5)
  return(p)
}
ggplot_peak_image <- function(pks,i,title="",col=1)
{
  palette=c("#8CD3C7","#FDFFB2","#BFBAD9")
  df=data.frame(x=pks$pos[,2],y=pks$pos[,1],z=pks$intensity[,i])
  background=element_rect(fill = palette[col])
  p=ggplot(df, aes(x, y, fill = z)) + geom_raster() +
    coord_fixed(1,expand = F) +
    ggtitle(title) +
    scale_fill_gradientn(colours=viridis(1000))+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_rect(fill = palette[col]),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=background,plot.title = element_text(hjust = 0.5,size=9))
  return(p)
}

plot_manual_validation<-function(x){
  res<-x$res
  table<-x$table
  plot(density(subset(res,adduct_mode=="neg"&fragmentation_possible&match)$lipid_occurences))
  lines(density(subset(res,adduct_mode=="neg"&fragmentation_possible&!match)$lipid_occurences),col=2)
  plot(density(subset(res,adduct_mode=="neg"&fragmentation_possible&match)$correlation))
  lines(density(subset(res,adduct_mode=="neg"&fragmentation_possible&!match)$correlation),col=2)
  plot(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&match),ambiguity*lipid_occurences)))
  lines(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&!match),ambiguity*lipid_occurences)),col=2)
  plot(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&match),correlation*lipid_occurences)))
  lines(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&!match),correlation*lipid_occurences)),col=2)
}
plot_manual_percentage <- function(res,scores,steps=100)
{
  table=load_gt()

  for(s in scores)
  {
    o<-rep(0,steps)
    th<-seq(0,max(res[,s]),length.out=steps)
    for(i in seq_along(th)){
      m<-manual_validation(res[res[,s]>th[i],],table)
      o[i]<-sum(m$table$match)/nrow(m$table)
    }
    plot(th,o)
  }
  return(o)
}

plot_manual_percentage_2 <- function(res,scores,steps=100)
{
  table=load_gt()

  for(s in scores)
  {
    o<-rep(0,steps)
    th<-seq(0,nrow(res),length.out=steps)
    a<-order(res[,s])
    for(i in seq_along(th)){
      m<-manual_validation(res[a[1:th[i]],],table)
      o[i]<-sum(m$table$match)/nrow(m$table)
    }
    plot(th,o)
  }
  return(o)
}
