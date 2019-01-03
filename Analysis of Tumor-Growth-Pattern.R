
##Copyright to D. H. Heiland

require(RColorBrewer)

# open the folder with files an example is given in the GITHUB

fileset=dir()

Set_Imaging=setClass("Set_Imaging",slots=c(
  Patdata="data.frame", 
  files="character",
  Time_dep_Frame="data.frame", 
  Growth="data.frame",
  Time_dep_Frame_Norm="data.frame",
  Data_Prop="data.frame"))


setMethod("initialize",signature = "Set_Imaging", definition = function(.Object, fileset ){
  .Object@files <- fileset
  return(.Object)
})

Set=Set_Imaging(fileset)



# Start by combining all patients into a time dependend matrix 

#Initiate virtual data.frame

data_time=data.frame(seq(0,7000, by=1)); names(data_time)="Time"
Set@Time_dep_Frame=data_time
Set@Time_dep_Frame_Norm=data_time
Set@Growth=data_time

#Import ime Frames into the virtual space
col=rainbow(length(Set@files))

i=1
for(i in i:length(Set@files)){
  print(i)
  in_dat=read.csv(Set@files[i], sep=";", header=T, dec=",")
  in_dat=in_dat[is.na(in_dat$Tumorvolumen)==F, ]
  Set@Time_dep_Frame[,i+1]=NA
  Set@Time_dep_Frame[Set@Time_dep_Frame$Time %in% in_dat$Zeit,i+1]=in_dat$Tumorvolumen
  if(i==1){
    plot(y=na.omit(Set@Time_dep_Frame[,i+1]),x=Set@Time_dep_Frame[is.na(Set@Time_dep_Frame[,i+1])==F,1],
         type="l", bty="n", xlim=c(0,7000), ylim=c(-10,150),
         main="Time Plot All Pat",
         xlab="Time in Days",
         ylab="Tumorvolume in ml")
    
  } else{
    points(y=na.omit(Set@Time_dep_Frame[,i+1]),x=Set@Time_dep_Frame[is.na(Set@Time_dep_Frame[,i+1])==F,1],type="l", col=col[i])
    
  }
  names(Set@Time_dep_Frame)[i+1]=Set@files[i]
  
  
}

#Check if all Pat in virtual Space
dim(Set@Time_dep_Frame)


#Normalize the size by start volume

i=1
for(i in i:length(Set@files)){
  print(i)
  in_dat=read.csv(Set@files[i], sep=";", header=T, dec=",")
  in_dat=in_dat[is.na(in_dat$Tumorvolumen)==F, ]
  in_dat$Tumorvolumen=in_dat[1,]$Tumorvolumen/in_dat$Tumorvolumen
  Set@Time_dep_Frame_Norm[,i+1]=NA
  Set@Time_dep_Frame_Norm[Set@Time_dep_Frame_Norm$Time %in% in_dat$Zeit,i+1]=in_dat$Tumorvolumen
  if(i==1){
    plot(y=na.omit(Set@Time_dep_Frame_Norm[,i+1]),x=Set@Time_dep_Frame_Norm[is.na(Set@Time_dep_Frame_Norm[,i+1])==F,1],
         type="l", bty="n", xlim=c(0,7000), ylim=c(-1,4),
         main="Time Plot All Pat",
         xlab="Time in Days",
         ylab="Fold Change in Tumor Volume")
    
  } else{
    points(y=na.omit(Set@Time_dep_Frame_Norm[,i+1]),x=Set@Time_dep_Frame_Norm[is.na(Set@Time_dep_Frame_Norm[,i+1])==F,1],type="l", col=col[i])
    
  }
  names(Set@Time_dep_Frame_Norm)[i+1]=Set@files[i]
  
  
  
  
  
  
  
  
}



#Validate tumor growth 
i=1
for(i in i:length(Set@files)){
  
  print(i)
  in_dat=read.csv(Set@files[i], sep=";", header=T, dec=",")
  in_dat=in_dat[is.na(in_dat$Tumorvolumen)==F, ]
  growth=na.omit(as.numeric(unlist(lapply(1:nrow(in_dat), function(i){
    growth=in_dat[i,]$Tumorvolumen/in_dat[i+1, ]$Tumorvolumen
    time_cor=in_dat[i+1, ]$Zeit/365
    cor_growth=growth/time_cor*100
  }))))
  
  in_dat$growth=c(0,as.numeric(growth))
  Set@Growth[,i+1]=NA
  Set@Growth[Set@Growth$Time %in% in_dat$Zeit,i+1]=in_dat$growth
  
  
  #plot growth
  
  if(i==1){
    plot(y=na.omit(Set@Growth[,i+1]),x=Set@Growth[is.na(Set@Growth[,i+1])==F,1],
         type="l", bty="n", xlim=c(0,7000), ylim=c(-1,1000),
         main="Time Plot All Pat",
         xlab="Time in Days",
         ylab="Fold Change in Tumor Volume")
    
  } else{
    points(y=na.omit(Set@Growth[,i+1]),x=Set@Growth[is.na(Set@Growth[,i+1])==F,1],type="l", col=col[i])
    
  }
  
  
  plot(y=na.omit(Set@Growth[,i+1]+20),x=Set@Growth[is.na(Set@Growth[,i+1])==F,1],
       type="l", bty="n", xlim=c(0,7000), ylim=c(-1,100),
       main="Time Plot All Pat",
       xlab="Time in Days",
       ylab="Fold Change in Tumor Volume")
  points(y=in_dat$Tumorvolumen,x=Set@Growth[is.na(Set@Growth[,i+1])==F,1],
         type="l", col="red", lwd=3)
  if (length(na.omit(in_dat$VolumenKM))>0) {points(y=in_dat$VolumenKM*10+20,x=Set@Growth[is.na(Set@Growth[,i+1])==F,1],type="l", col="green", lwd=3)}
  
  
  
  
  
  require(e1071)
  
  #Validate Tumor growth
  
  dens=density(growth)
  min_grow=min(growth)
  max_grow=max(growth)
  medium_grow=mean(growth)
  median_grow=median(growth)
  variation_grow=var(growth)
  kurtosis_grow=kurtosis(growth)
  skrewness_grow=skewness(growth)
  dens=density(in_dat$Tumorvolumen)
  min_V=min(in_dat$Tumorvolumen)
  max_V=max(in_dat$Tumorvolumen)
  medium_V=mean(in_dat$Tumorvolumen)
  median_V=median(in_dat$Tumorvolumen)
  variation_V=var(in_dat$Tumorvolumen)
  kurtosis_V=kurtosis(in_dat$Tumorvolumen)
  skrewness_V=skewness(in_dat$Tumorvolumen)
  
  grow_prop=c(min_grow,max_grow,medium_grow,median_grow,variation_grow,kurtosis_grow,skrewness_grow,min_V,max_V,medium_V,median_V,variation_V,kurtosis_V,skrewness_V)
  names_grow=c("min_grow","max_grow","medium_grow","median_grow","variation_grow","kurtosis_grow","skrewness_grow","min_V","max_V","medium_V","median_V","variation_V","kurtosis_V","skrewness_V")
  
  out=data.frame(grow_prop)
  names(out)=paste(names(out),Set@files[i], sep="_")
  rownames(out)=names_grow
  names(Set@Growth)[i+1]=Set@files[i]
  
  
  ### Add response to data
  lv=c("PCV","PC","Temodal","CCNU", "Resektion", "RTX")
  
  for (ix in 1: length(lv)){
    #test if Drug is avaiable
    print(ix)
    test=length(na.omit(in_dat[in_dat[,lv[ix]]==1,lv[ix]]))
    
    
    
    
    
    if (test>=2){
      rownames(in_dat)=seq(1:nrow(in_dat))
      db=in_dat[in_dat[,lv[ix]]==1,]
      times=rownames(db[is.na(db[,lv[ix]])==F, ])
      if(as.numeric(times[length(times)])==length(rownames(db))){times=c(as.numeric(times[length(times)])-1, times)}else{ times=c(times,as.numeric(times[length(times)])+1)}
      
      
      Tumorvolumen=in_dat[times, ]$Tumorvolumen
      growth=in_dat[times, ]$growth
      
      if(lv[ix]=="PC"){Response_PCV[i,]=as.numeric(Tumorvolumen[1:3])}
      if(lv[ix]=="Resektion"){Response_Re[i, ]=as.numeric(Tumorvolumen[1:3])}
      if(lv[ix]=="RTX"){Response_RTX[i, ]=as.numeric(Tumorvolumen[1:3])}
      
      
      
      
      
      
      #Res_dens=density(growth)
      Res_min_grow=min(growth)
      Res_max_grow=max(growth)
      Res_medium_grow=mean(growth)
      Res_median_grow=median(growth)
      Res_variation_grow=var(growth)
      Res_kurtosis_grow=kurtosis(growth)
      Res_skrewness_grow=skewness(growth)
      #Res_dens=density(Tumorvolumen)
      Res_min_V=min(Tumorvolumen)
      Res_max_V=max(Tumorvolumen)
      Res_medium_V=mean(Tumorvolumen)
      Res_median_V=median(Tumorvolumen)
      Res_variation_V=var(Tumorvolumen)
      Res_kurtosis_V=kurtosis(Tumorvolumen)
      Res_skrewness_V=skewness(Tumorvolumen)
      
      
      
      grow_prop_x=c(Res_min_grow,Res_max_grow,
                    Res_medium_grow,Res_median_grow,
                    Res_variation_grow,Res_kurtosis_grow,
                    Res_skrewness_grow,Res_min_V,Res_max_V,
                    Res_medium_V,Res_median_V,Res_variation_V,
                    Res_kurtosis_V,Res_skrewness_V)
      names_grow_x=c("Res_min_grow","Res_max_grow","Res_medium_grow","Res_median_grow","Res_variation_grow",
                     "Res_kurtosis_grow","Res_skrewness_grow","Res_min_V","Res_max_V","Res_medium_V",
                     "Res_median_V","Res_variation_V","Res_kurtosis_V","Res_skrewness_V")
      
      
      
      
      out_Resp=data.frame(grow_prop_x)
      rownames(out_Resp)=paste(names_grow_x, lv[ix], sep="_")
      names(out_Resp)=names(out)
      
      out=rbind(out,out_Resp)
      
      
      
      
    }
    if (test<2){
      out_Resp=data.frame(rep(0, 14))
      rownames(out_Resp)=paste(names_grow_x, lv[ix], sep="_")
      names(out_Resp)=names(out)
      
      out=rbind(out,out_Resp)
      
    }
    
    
    
    
    
  }
  
  
  
  
  
  print(dim(out))
  
  if(i==1){Set@Data_Prop=out}else{Set@Data_Prop=cbind(Set@Data_Prop,out)}
  
  
  
  
  
  
  
}
Set@Data_Prop

y=(scale(as.matrix(Set@Data_Prop)))


#potential plots for response

Response_PCV=data.frame(row.names=PIZ, R1=rep(1,44),R2=rep(1,44),R3=rep(1,44))
Response_RTX=data.frame(row.names=PIZ, R1=rep(1,44),R2=rep(1,44),R3=rep(1,44))
Response_Re=data.frame(row.names=PIZ, R1=rep(1,44),R2=rep(1,44),R3=rep(1,44))

data_res=Response_PCV

data_res=data_res[!rowSums(data_res[, 1:3 ])==3, ]
for(i in 1:nrow(data_res)){data_res[i, ]=data_res[i,]/data_res[i,1]}
data_res$col="green"
data_res[c(data_res[,2]-data_res[,1 ])>0, ]$col="red"

plot(x=c(1,2,3), y=data_res[1,1:3], type="l", ylim=c(0,2), bty="n", col=data_res$col[1])
for(i in 2:nrow(data_res)){points(x=c(1,2,3), y=data_res[i,1:3], type="l", col=data_res$col[i])}




#hier eine Heatmap mit features erstellen

clinic_data=read.csv("clinical.csv", row.names = 1, sep=";")
Set@Patdata=clinic_data
PIZ=NULL
for (i in 1:length(Set@files)) {a=as.numeric(unlist(strsplit(unlist(strsplit( colnames(y)[i], "_"))[4],".csv"))[1]);print(a);PIZ=c(PIZ,a) }

groups_men=groups_men[names(Set@Data_Prop), ]

Set@Patdata=Set@Patdata[as.character(PIZ), ]
Set@Patdata$col_WHO="red"

#.... add more feratures....


heatmap.plus::heatmap.plus(cor(y), col=brewer.pal(9,"RdPu" ),labRow=T,labCol=T,ColSideColors=as.matrix(cbind(
  Set@Patdata$col_WHO
  #... add more
)))

heatmap(cor((y)), col=rev(brewer.pal(9,"Spectral" )))

#Ausarbeiten welches feature besonders mit welcher gruppe assoziiert ist
res<-AutoPipe::TopPAM((y),max_clusters = 8, TOP=98)
me_TOP=res[[1]]


