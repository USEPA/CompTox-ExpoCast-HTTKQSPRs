                           # This function does the level II concentration comparisons:
makeCvTpreds <- function(label)
{
  cvt.table <- NULL
  stats.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
  {
    if (this.cas %in% get_cheminfo(model="pbtk"))
    {
      this.subset1 <- subset(CvT.data,CAS==this.cas & 
        !is.na(Time) & 
        !is.na(Value))
      this.compound <- this.subset1$Compound[1]
      this.dtxsid <- this.subset1$DTXSID[1]
      for (this.species in unique(this.subset1$Species))
      {
        this.subset2 <- subset(this.subset1,Species==this.species)
        for (this.route in unique(this.subset2$Route))
        {
          this.subset3 <- subset(this.subset2,Route==this.route)
          for (this.dose in unique(this.subset3$Dose))
          {
            this.subset4 <- subset(this.subset3,Dose==this.dose)
            obs.times <- signif(sort(unique(this.subset4$Time)),4)
            pred <- suppressWarnings(solve_pbtk(
              chem.cas=this.cas,
              times=sort(unique(c(seq(0,2,0.05),obs.times))),
              species=this.species,
              iv.dose=(this.route=="iv"),
              dose=this.dose,
              default.to.human=TRUE,
              suppress.messages=TRUE))
            pred <- subset(pred,pred[,"time"]>0.0001)
            # Convert from uM to ug/mL:
            pred[,"Cven"] <- pred[,"Cven"]*parameterize_pbtk(chem.cas=this.cas)$MW/1000
            if (any(tolower(unlist(this.subset4[,"Media"]))=="blood"))
            {
              Rb2p <- available_rblood2plasma(
                chem.cas=this.cas,
                species=this.species)
            }
            this.subset4means <- NULL
            for (this.time in obs.times)
            {
              this.subset5 <- subset(this.subset4,signif(Time,4)==this.time)
  #            print(this.subset5)
              new.row <- data.frame(
                Compound=this.compound,
                DTXSID=this.dtxsid,
                CAS=this.cas,
                Species=this.species,
                Route=this.route,
                Dose=this.dose,
                Time=this.time,
                Conc.pred=pred[pred[,"time"] == this.time,"Cven"],
                stringsAsFactors=F)
              new.row <- merge(new.row,
                this.subset5[,c(
                  "CAS","Media","Value","Source","calc_loq")], by="CAS")
              colnames(new.row)[colnames(new.row)=="Value"] <- "Conc.obs"
              if (any(tolower(new.row[,"Media"])=="blood"))
              {
                new.row[tolower(new.row[,"Media"])=="blood","Conc.pred"] <-
                  Rb2p * new.row[tolower(new.row[,"Media"])=="blood","Conc.pred"]
              }
              cvt.table <- rbind(cvt.table,new.row)
              this.subset4means <- rbind(this.subset4means,
                data.frame(
                  Time=this.time,
                  Value=mean(this.subset5$Value,na.rm=0)))
            }
            Cmax.obs <- max(this.subset4means$Value,na.rm=T)
            Cmax.pred <- max(pred[,"Cven"],na.rm=T)
            AUC.obs <- AUC(this.subset4means$Time,this.subset4means$Value)
            AUC.pred <- AUC(pred[,"time"],pred[,"Cven"])
            new.row <- data.frame(
              Compound=this.compound,
              DTXSID=this.dtxsid,
              CAS=this.cas,
              Species=this.species,
              Route=this.route,
              Dose=this.dose,
              Cmax.obs=Cmax.obs,
              Cmax.pred=Cmax.pred,
              AUC.obs=AUC.obs,
              AUC.pred=AUC.pred,
              stringsAsFactors=F)
            stats.table <- rbind(stats.table,new.row)            
          }
        }
      }
    }
  }
  cvt.table$QSAR <- label
  stats.table$QSAR <- label
  return(list(
    cvt=cvt.table,
    stats=stats.table))
}


# This function does the level II concentration comparisons:
makeCvTpredsfromfits <- function(label)
{
  cvt.table <- NULL
  stats.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
  {
    if (this.cas %in% cvtfits$CAS & 
      this.cas %in% get_cheminfo(model="1compartment"))
    {
      this.fit <- subset(cvtfits, CAS==this.cas)
      this.subset1 <- subset(CvT.data,CAS==this.cas & 
        !is.na(Time) & 
        !is.na(Value))
      this.compound <- this.subset1$Compound[1]
      this.dtxsid <- this.subset1$DTXSID[1]
      for (this.species in unique(this.subset1$Species))
        if (tolower(this.species) %in% tolower(this.fit$Species))
        {
          params <- suppressWarnings(parameterize_1comp(
            chem.cas=this.cas,
            species=this.species,
            default.to.human=TRUE,
            suppress.messages=TRUE))
          fit.index <- tolower(this.fit$Species)==tolower(this.species)
          params$Vdist <- as.numeric(this.fit[fit.index,"Vdist"])
          params$kelim <- as.numeric(this.fit[fit.index,"kelim"]) 
          if (!is.na(this.fit[fit.index,"Fgutabs"])) 
          {
            params$Fgutabs <- as.numeric(this.fit[1,"Fgutabs"])
          }
          if (!is.na(this.fit[fit.index,"kgutabs"])) 
          {
            params$kgutabs<- as.numeric(this.fit[1,"kgutabs"])
          }
          this.subset2 <- subset(this.subset1,Species==this.species)
          for (this.route in unique(this.subset2$Route))
          {
            this.subset3 <- subset(this.subset2,Route==this.route)
            for (this.dose in unique(this.subset3$Dose))
            {
              this.subset4 <- subset(this.subset3,Dose==this.dose)
              obs.times <- signif(sort(unique(this.subset4$Time)),4)
              pred <- suppressWarnings(solve_1comp(
                parameters=params,
                times=sort(unique(c(seq(0,2,0.05),obs.times))),
                iv.dose=(this.route=="iv"),
                dose=this.dose,
                suppress.messages=TRUE))
              pred <- subset(pred,pred[,"time"]>0.0001)
              # Convert from uM to ug/mL:
              pred[,"Ccompartment"] <- pred[,"Ccompartment"]*params$MW/1000
              if (any(tolower(unlist(this.subset4[,"Media"]))=="blood"))
              {
                Rb2p <- suppressWarnings(available_rblood2plasma(
                  chem.cas=this.cas,
                  species=this.species,
                  suppress.messages=TRUE))
              }
              this.subset4means <- NULL
              for (this.time in obs.times)
              {
                this.subset5 <- subset(this.subset4,signif(Time,4)==this.time)
                new.row <- data.frame(
                  Compound=this.compound,
                  DTXSID=this.dtxsid,
                  CAS=this.cas,
                  Species=this.species,
                  Route=this.route,
                  Dose=this.dose,
                  Time=this.time,
                  Conc.pred=pred[pred[,"time"] == this.time,"Ccompartment"],
                  stringsAsFactors=F)
                new.row <- merge(new.row,
                  this.subset5[,c(
                    "CAS","Media","Value","Source","calc_loq")], by="CAS")
                colnames(new.row)[colnames(new.row)=="Value"] <- "Conc.obs"
                if (any(tolower(new.row[,"Media"])=="blood"))
                {
                  new.row[tolower(new.row[,"Media"])=="blood","Conc.pred"] <-
                    Rb2p * new.row[tolower(new.row[,"Media"])=="blood","Conc.pred"]
                }
                cvt.table <- rbind(cvt.table,new.row)
                this.subset4means <- rbind(this.subset4means,
                  data.frame(
                    Time=this.time,
                    Value=mean(this.subset5$Value,na.rm=0)))
              }
              Cmax.obs <- max(this.subset4means$Value,na.rm=T)
              Cmax.pred <- max(pred[,"Ccompartment"],na.rm=T)
              AUC.obs <- AUC(this.subset4means$Time,this.subset4means$Value)
              AUC.pred <- AUC(pred[,"time"],pred[,"Ccompartment"])
              new.row <- data.frame(
                Compound=this.compound,
                DTXSID=this.dtxsid,
                CAS=this.cas,
                Species=this.species,
                Route=this.route,
                Dose=this.dose,
                Cmax.obs=Cmax.obs,
                Cmax.pred=Cmax.pred,
                AUC.obs=AUC.obs,
                AUC.pred=AUC.pred,
                stringsAsFactors=F)
              stats.table <- rbind(stats.table,new.row)            
            }
          }
      }
    }
  }
  cvt.table$QSAR <- label
  stats.table$QSAR <- label
  return(list(
    cvt=cvt.table,
    stats=stats.table))
}

# This function does the level III statistic comparisons:
maketkstatpreds <- function(label)
{
  out.table <- NULL
  for (this.cas in unique(cvtfits$CAS))
  {
    if (this.cas %in% get_cheminfo(model="1compartment"))
    {
      this.subset1 <- subset(cvtfits,CAS==this.cas)
      this.compound <- this.subset1$Compound[1]
      this.dtxsid <- this.subset1$DTXSID[1]
      for (this.species in unique(this.subset1$Species))
      {
        this.subset2 <- subset(this.subset1,Species==this.species)
        vd.obs <- unlist(as.numeric(this.subset2[,"Vdist"]))
        vd.pred <- suppressWarnings(calc_vdist(
          chem.cas=this.cas,
          species=this.species,
          default.to.human=TRUE,
          suppress.messages=TRUE))
        thalf.obs <- unlist(as.numeric(this.subset2[,"halflife"]))
        thalf.pred <- log(2)/suppressWarnings(calc_elimination_rate(
          chem.cas=this.cas,
          species=this.species,
          default.to.human=TRUE,
          suppress.messages=TRUE))
        new.tab <- data.frame(
          Compound=this.compound,
          DTXSID=this.dtxsid,
          CAS=this.cas,
          Species=this.species,
          Vd.obs=vd.obs,
          Vd.pred = vd.pred,
          thalf.obs=thalf.obs,
          thalf.pred = thalf.pred,
          stringsAsFactors=F)
        new.tab <- merge(new.tab,this.subset2[,
          c("CAS","Reference")],by="CAS")
        out.table <- rbind(out.table,new.tab)
      }
    }
  }
  out.table$QSAR <- label
  return(out.table)
}
