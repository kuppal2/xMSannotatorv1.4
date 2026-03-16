
library(xMSannotator)

###########Parameters to change##############
system.time(annotres4<-multilevelannotation(
  max.mz.diff=10, max.rt.diff=10,corthresh=0.7,
                        dataA=read.csv(system.file("extdata","StdMix_feature_table.csv",package = "xMSannotator")),
                       outloc="/Users/karanuppal//xMSannotatorv1.4/", 
                          mode="pos",
                                            clustmethod="auto", 
                                            queryadductlist=c("M+H","M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O"),
                                            db_name="HMDB",
                                            customIDs=read.csv(system.file("extdata","StdMix_HMDBIDs.csv",package = "xMSannotator")),
                                            missing.value=0,
                                            filter.by=c("M+H"),
                                             boostIDs=NA,
                                             max_isp=5,
                                            customDB=NA,
                                            peakID_name = "PeakID"
  )
)

