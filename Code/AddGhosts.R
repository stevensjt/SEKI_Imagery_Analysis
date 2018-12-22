AddGhosts <- function(SoilM){
  
  AllYrs<-unique(SoilM$Year)
  AllSubSites<-sort(unique(SoilM$SubSite))
  
  Count1e=0;
  Count1l=0;
  Count2e=0;
  Count2l=0;
  
  for(i in AllSubSites){
    i
    npts<-length(SoilM$VWC)
    SubMat<-SoilM[SoilM$SubSite==i,]
    YrsEarlySummer<-unique(SubMat$Year[SubMat$DOY<180])
    NYES<-length(YrsEarlySummer)
    YrsLateSummer<-unique(SubMat$Year[SubMat$DOY>180])
    NYLS<-length(YrsLateSummer)
    #if(NYLS==0){}
    #if(NYES==0){}
    #How should I fill in data? Linear rule?
    #Should I make the cutoff 0 or 1?
    if(NYES!=NYLS){
      if ((NYLS==(NYES-1))&&(NYLS==1)) {#Replicate year of late summer data
        
        Mssng<-which(YrsEarlySummer%in%YrsLateSummer)
        if (length(Mssng)>0){
          Temp1=SubMat[(SubMat$DOY>180),]
          SoilM[(npts+1):(npts+length(Temp1$DOY)),]=Temp1
          SoilM$Year[(npts+1):(npts+length(Temp1$DOY))]=YrsEarlySummer[-Mssng]
          Count1l=Count1l+1
          }
        else{#NOT DONE. What if there are 2 years w/ both early and late summer but one w/out? Is it worth filling it in?
        }
        npts=length(SoilM$VWC)
      }
      else if ((NYES==(NYLS-1))&&(NYES==1)) {#Replicate year of early summer data
        Mssng<-which(YrsLateSummer%in%YrsEarlySummer)
        if (length(Mssng)>0){
          Count1e=Count1e+1
          Temp1=SubMat[(SubMat$DOY<180),]
          SoilM[(npts+1):(npts+length(Temp1$DOY)),]=Temp1
          SoilM$Year[(npts+1):(npts+length(Temp1$DOY))]=YrsLateSummer[-Mssng]
        }
        else{#NOT DONE
        }
        npts=length(SoilM$VWC)
      }
      else if (NYES==0) {
        
        if (mean(SubMat$VWC)<20){mvwc=2.9}else{mvwc=1.4}
        NYES
        mvwc
        for (j in YrsLateSummer) {
          Count2e=Count2e+1
          LateSumm1=SubMat[(SubMat$Year==j),]
          SoilM[(npts+1):(npts+length(LateSumm1$DOY)),]=LateSumm1
          SoilM$Year[(npts+1):(npts+length(LateSumm1$DOY))]=j
          SoilM$DOY[(npts+1):(npts+length(LateSumm1$DOY))]=150
          SoilM$VWC[(npts+1):(npts+length(LateSumm1$DOY))]=pmin(55,mvwc*LateSumm1$VWC)
          npts=length(SoilM$VWC)              
        }}
      else if (NYLS==0) {
        
        if (mean(SubMat$VWC)<30){mvwc=.35
         #If it's a relatively dry site
        } else {mvwc=.8 #wetter sites were between .6 and 1 of early summer VWC in late summer
        }
        NYLS
        mvwc
        for (j in YrsEarlySummer) {
          Count2l=Count2l+1
          EarlySumm1=SubMat[(SubMat$Year==j),]
          SoilM[(npts+1):(npts+length(EarlySumm1$DOY)),]=EarlySumm1
          SoilM$Year[(npts+1):(npts+length(EarlySumm1$DOY))]=j
          SoilM$DOY[(npts+1):(npts+length(EarlySumm1$DOY))]=210
          SoilM$VWC[(npts+1):(npts+length(EarlySumm1$DOY))]=mvwc*EarlySumm1$VWC #drier sites decreased at about 0.001/day or 1-2%/day. Mean proportion of late to early was around 0.35.
          npts=length(SoilM$VWC)
        }
      }
    }
  }
  
  return(SoilM)
}