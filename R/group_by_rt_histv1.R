group_by_rt_histv1 <-
function(mchemicaldata,time_step=2,max_diff_rt=10,groupnum){
  
  
  range1<-range(mchemicaldata$time)
  diff1<-abs(range1[1]-range1[2])
  
  if(isTRUE((nrow(mchemicaldata)>1 && diff1>1)[1])){
    
    d1<-try(density(mchemicaldata$time,bw=time_step)) #,adjust=0.1))
    ts_y<-ts(d1$y)
    
    tp<-getturnpoints(ts_y)
    #find.peak.regions(ts_y,nups=1,threshold=(-1),minpeakheight=(-1),minpeakdistance=time_step) #
    
    pdfname<-paste("timeclust_groupnum.pdf",sep="")
    pdf(pdfname)
    plot(d1)
    points(d1$x[tp$tppos],d1$y[tp$tppos],col="red")
    dev.off()
    
    #points(d1$x[tp$tppos][seq(1,length(tp$tppos),2)],d1$y[tp$tppos][seq(1,length(tp$tppos),2)],col="red")
    
    
    time_break_points<-d1$x[tp$tppos] #[seq(1,length(tp$tppos),2)]
    
    if (isTRUE((is(d1, "try-error"))[1])){
      bw_i<-max_diff_rt
    }else{
      
      
      bw_i<-d1$bw
    }
    
  }else{
    bw_i<-max_diff_rt
    time_break_points<-mchemicaldata$time
    
  }
  
  
  
  
  
  t1<-mchemicaldata$time
  t1<-as.data.frame(t1)
  colnames(t1)<-c("mids")
  
  #time_cor_groups<-sapply(list(myData1=bh$x),function(x)  split(x,cut(bh$x,breaks=bh$x[tp$tppos])))
  
  
  # time_cor_groups<-sapply(list(myData1=h1$density),function(x)  split(x,cut(h1$mids,breaks=seq(min(h1$mids)-max_diff_rt,max(h1$mids)+max_diff_rt,10))))
  clusternum<-1;clusterlabels=rep(1,dim(t1)[1]);
  
  for(i in seq(1,length(time_break_points),2)){
    
    
    if(isTRUE(((i+1)>length(time_break_points))[1])){
      
      clusterlabels[which(t1$mids>=time_break_points[i-1])]<-clusternum
    }else{
      
      if(isTRUE((i==1)[1])){
        
        clusterlabels[which(t1$mids<=time_break_points[i] & t1$mids>=time_break_points[i] & t1$mids<time_break_points[i+1])]<-clusternum
      }else{
        clusterlabels[which(t1$mids>=time_break_points[i-1] & t1$mids<time_break_points[i+1])]<-clusternum
      }
    }
    
    clusternum<-clusternum+1
    
  }
  t2<-cbind(t1,clusterlabels)
  t2<-as.data.frame(t2)
  
  
  
  levelBnum<-paste(groupnum,clusterlabels,sep="_")
  
  diffmatB<-cbind(levelBnum,mchemicaldata)
  
  
  
  diffmatB<-as.data.frame(diffmatB)
  
  
  cnames<-colnames(diffmatB)
  cnames[1]<-"Module_RTclust"
  colnames(diffmatB)<-cnames
  
  return(diffmatB)
  
  
}
