                           # This function does the level II concentration comparisons:
makeCvTpreds <- function(CvT.data,label)
{
  nonvol.chems <- get_cheminfo(model="pbtk")
  vol.chems <- get_cheminfo(model="gas_pbtk")
  
  cvt.table <- NULL
  stats.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
  {
    this.model <- NULL
    if (this.cas %in% nonvol.chems) this.model <- "solve_pbtk"
    else if (this.cas %in% vol.chems) this.model <- "solve_gas_pbtk"
    
    if (!is.null(this.model))
    {
 #     print(this.cas)
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

            pred <- suppressWarnings(eval(call(this.model,
              chem.cas=this.cas,
              times=sort(unique(c(seq(0,2,0.05),obs.times))),
              species=this.species,
              iv.dose=(this.route=="iv"),
              dose=this.dose,
              default.to.human=TRUE,
              suppress.messages=TRUE,
              input.units='mg/kg',
              exp.conc=0)))

            pred <- subset(pred,pred[,"time"]>0.0001)
            # Convert from uM to ug/mL:
            pred[,"Cven"] <- pred[,"Cven"]*suppressWarnings(parameterize_pbtk(
              chem.cas=this.cas,
              suppress.messages=TRUE))$MW/1000
            if (any(tolower(unlist(this.subset4[,"Media"]))=="blood"))
            {
              Rb2p <- suppressWarnings(available_rblood2plasma(
                chem.cas=this.cas,
                species=this.species,
                suppress.messages = TRUE))
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
            if (length(unique(this.subset4means$Time))>1)
            {
              AUC.obs <- AUC(this.subset4means$Time,this.subset4means$Value)
              AUC.pred <- AUC(pred[,"time"],pred[,"Cven"])
            } else {
              AUC.obs <- NA
              AUC.pred <- NA
            }
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
  cvt.table$QSPR <- label
  stats.table$QSPR <- label
  return(list(
    cvt=cvt.table,
    stats=stats.table))
}

#'Analytical 1-compartment model
#'
#'@param time A vector of times in hours
#'@param params A named list of model parameter values. Must include:
#'\describe{
#'\item{kelim}{Elimination rate, 1/h}
#'\item{Vdist}{Volume of distribution, L/kg body weight}}
#'For oral administration (\code{iv.dose} FALSE), \code{params} must also include:
#'\describe{
#'\item{Fgutabs}{Oral bioavailability, unitless fraction}
#'\item{kgutabs}{Oral absorption rate, 1/h}}
#'@param dose Dose in mg/kg
#'@param iv.dose TRUE for single IV bolus dose; FALSE for single oral dose
#'
#'@return A vector of plasma concentration values corresponding to \code{time}.
#'
#'@author Caroline Ring, John Wambaugh
#'
#' @export cp_1comp
cp_1comp <- function(time, params, dose, iv.dose){
  #time in hours
  #dose in mg/kg
  #params: a subset of those returned by httk::parameterize_1comp
  #a named list
  #Fgutabs = oral fraction absorbed, unitless
  #kelim = elimination rate, 1/h
  #kgutabs = oral absorption rate, 1/h
  #Vdist = volume of distribution, L/kg body weight
  if (any(sapply(params,function(x) identical(x,numeric(0))))) return(0)
  
  if (is.null(params$Fgutabs)|is.na(params$Fgutabs)) 
  {
    params$Fgutabs<-1
  }
  if (is.null(params$kgutabs)|is.na(params$kgutabs)) 
  {
    params$kgutabs<-1
  }  
  if (params$Fgutabs>1) params$Fgutabs<-1
  
  if (iv.dose){
    cp <- dose*exp(-params$kelim * time)/params$Vdist
  }else{
    
    cp <- (params$Fgutabs * dose *
             params$kgutabs*(exp(-params$kelim * time) - exp(-params$kgutabs* time)))/
      (params$Vdist*(params$kgutabs - params$kelim))
    #Note: this fails if kgutabs == kelim,
    #which usually happens if both are at the upper bound
  }
  cp[cp<0] <- 0
  cp[!is.finite(cp)] <- 0
  return(cp)
}

#' Analytical 2-compartment model
#' 
#' 
#' @param params A named list of parameter values including the following:
#' \describe{ \item{k12}{Rate at which compound moves from central to
#' peripheral compartment} \item{k21}{Rate at which compound moves from
#' peripheral to central compartment} \item{kelim}{Elimination rate}
#' \item{V1}{Apparent volume of central compartment} } For oral administration
#' (\code{iv.dose} FALSE), \code{params} must also include: \describe{
#' \item{Fgutabs}{Oral bioavailability} \item{kgutabs}{Rate of absorption from
#' gut} }
#' 
#' !author Caroline Ring, John Wambaugh
#' @param time A vector of time values, in hours
#' @param dose A dose in mg/kg
#' @param iv.dose TRUE for single IV bolus dose, FALSE for single oral dose
#' @return A vector of plasma concentration values corresponding to each value
#' in \code{time}
#' @export cp_2comp
cp_2comp <- function(params, time, dose, iv.dose)
{
  if (any(sapply(params,function(x) identical(x,numeric(0))))) return(NA)
  
  if (is.null(params$Fgutabs)|is.na(params$Fgutabs)) 
  {
    params$Fgutabs<-1
  }
  if (is.null(params$kgutabs)|is.na(params$kgutabs)) 
  {
    params$kgutabs<-1
  }
  # if (is.null(params$Fbetaofalpha)) params$Fbetaofalpha<-1
  # if (is.null(params$Ralphatokelim)) params$Ralphatokelim<-1
  # 
  # if (is.na(params$Fbetaofalpha)) params$Fbetaofalpha<-1
  # if (is.na(params$Ralphatokelim)) params$Ralphatokelim<-1
  # 
  # if (params$Fgutabs>1) params$Fgutabs<-1
  # if (params$Fbetaofalpha>1) params$Fbetaofalpha<-1
  # if (params$Ralphatokelim<0) params$Ralphatokelim<-1
  # 
  # alpha <- params$Ralphatokelim*(params$kelim+10^-6)
  # beta <- params$Fbetaofalpha*alpha
  # 
  # # try to keep k21 and k12 positive:
  # k21 <- max(min(alpha*beta/params$kelim,alpha+beta-params$kelim),0)
  # k12 <- alpha + beta - params$kelim - k21
  
  # Elimination rate from the body:
  kelim <- params$kelim
  # Transfer rate from deep tissue to primary:
  k21 <- params$k21
  # Transfer rate from primary to deep tissue:
  k12 <- params$k12
  # Absorption rate from gut:
  kgutabs <- params$kgutabs
  # Primary compartment volume:
  V1 <- params$V1
  # Fraction systemically bioavailable from oral dose:
  Fgutabs <- params$Fgutabs
  
  alphabeta.sum <- k12 + k21 + kelim
  alphabeta.prod <- k21*kelim
 
  beta <- 1/2*(alphabeta.sum - (alphabeta.sum^2-4*alphabeta.prod)^(1/2))
  alpha <- 1/2*(alphabeta.sum + (alphabeta.sum^2-4*alphabeta.prod)^(1/2))   
  
  if (iv.dose){ #for IV dosing
    A <- (dose*(alpha - k21))/(V1*(alpha-beta))
    B <- (dose*(k21-beta))/(V1*(alpha-beta))
    C <- 0
  }else{ #for oral dosing
    A <- dose*(Fgutabs/V1)*(kgutabs/(kgutabs-alpha))*((alpha - k21)/(alpha-beta))
    B <- dose*(Fgutabs/V1)*(kgutabs/(kgutabs-beta))*((k21-beta)/(alpha-beta))
    C <- -(A+B)
  }

  cp <-  A*exp(-alpha*time) + B*exp(-beta*time) + C*exp(-kgutabs*time)
  cp[cp < 10^-20] <- 10^-20
  return(cp)
}

# This function does the level II concentration comparisons:
makeCvTpredsfromfits <- function(
  CvT.data, # Concentration vs. time data
  fittable,
  label = "FitsToData") # How the predictions should be labeled
{
  cvt.table <- NULL
  stats.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
  {
    this.subset1 <- subset(CvT.data,CAS==this.cas & 
                             !is.na(Time) & 
                             !is.na(Value))
    this.compound <- this.subset1$Compound[1]
    this.dtxsid <- this.subset1$DTXSID[1]
    
    if (this.cas %in% fittable$CAS)
    {
      this.fit1 <- subset(fittable, CAS==this.cas)

      for (this.species in unique(tolower(this.fit1$Species)))
      {
        this.fit2 <- subset(this.fit1,tolower(Species)==this.species)

        model <- NULL
        if (this.fit2$Model %in% c("1Comp","2Comp"))
        {
          if (this.fit2$Model == "1Comp") model <- "cp_1comp"
          else model <- "cp_2comp"
          
          params <- NULL
          fit.index <- 1
          if (model == "cp_1comp")
          {
            params$Vdist <- as.numeric(this.fit2[fit.index, "Vdist"])
          } else {
            params$V1 <- as.numeric(this.fit2[fit.index, "V1"])
            params$k12 <- as.numeric(this.fit2[fit.index,"k12"])
            params$k21 <- as.numeric(this.fit2[fit.index,"k21"])
          }
          
          # Same parameter names regardless of model:  
          params$kelim <- as.numeric(this.fit2[fit.index,"kelim"]) 
          params$Fgutabs <- as.numeric(this.fit2[fit.index, "Fgutabs"])
          params$kgutabs<- as.numeric(this.fit2[fit.index, "kgutabs"])
          
          this.subset2 <- subset(this.subset1,tolower(Species)==this.species)
          for (this.route in unique(this.subset2$Route))
          {
            this.subset3 <- subset(this.subset2,Route==this.route)
            for (this.dose in unique(this.subset3$Dose))
            {
              this.subset4 <- subset(this.subset3,Dose==this.dose)
              obs.times <- signif(sort(unique(this.subset4$Time)),4)
              # Dose is in mg/kg. Vd is in L/kg
              pred <- suppressWarnings(signif(do.call(
                model,
                list(
                  time = obs.times,
                  params = params,
                  dose = this.dose,
                  iv.dose = tolower(this.route)=="iv")), 4))
              if (any(is.na(pred))) browser()
              # I think output is already in ug/mL (mg/L)
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
                  Conc.pred=pred[obs.times == this.time],
                  stringsAsFactors=FALSE)
                new.row <- merge(new.row,
                                 this.subset5[,c(
                                   "CAS","Media","Value","Source","calc_loq")], by="CAS")
                colnames(new.row)[colnames(new.row)=="Value"] <- "Conc.obs"
                cvt.table <- rbind(cvt.table,new.row)
                this.subset4means <- rbind(this.subset4means,
                                           data.frame(
                                             Time=this.time,
                                             Value=mean(this.subset5$Value,na.rm=0)))
              }
              Cmax.obs <- max(this.subset4means$Value,na.rm=T)
              Cmax.pred <- max(pred,na.rm=T)
            #  if (Cmax.pred == 0) browser()
              if (length(obs.times)>1)
              {
                AUC.obs <- AUC(this.subset4means$Time,this.subset4means$Value)
                AUC.pred <- AUC(obs.times,pred)
              } else {
                AUC.obs <- NA
                AUC.pred <- NA
              }
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
                stringsAsFactors=FALSE)
              stats.table <- rbind(stats.table,new.row)            
            }
          }
        }
      }
    }
  }
        
  cvt.table$QSPR <- label
  stats.table$QSPR <- label
  return(list(
    cvt=cvt.table,
    stats=stats.table))
}

# This function does the level II concentration comparisons:
makeCvTpredsfromfits2 <- function(
    CvT.data, # Concentration vs. time data
    fit.preds,
    label = "FitsToData") # How the predictions should be labeled
{
  cvt.table <- NULL
  stats.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
  {
    this.cvt.subset1 <- subset(CvT.data,
                               CAS==this.cas & 
                               !is.na(Time) & 
                               !is.na(Value))
    this.pred.subset1 <- subset(fit.preds, 
                                CAS==this.cas & 
                                !is.na(Time) & 
                                !is.na(Value))
    this.compound <- this.cvt.subset1$Compound[1]
    this.dtxsid <- this.cvt.subset1$DTXSID[1]
    
    for (this.species in unique(tolower(this.cvt.subset1$Species)))
      if (this.species %in% tolower(this.pred.subset1$Species))
      {
        this.pred.subset2 <- subset(this.pred.subset1,
                                    tolower(Species)==this.species)
        this.cvt.subset2 <- subset(this.cvt.subset1,
                                    tolower(Species)==this.species)

        for (this.route in unique(tolower(this.cvt.subset2$Route)))
          if (this.route %in% tolower(this.pred.subset2$Route))
          {
            this.cvt.subset3 <- subset(this.cvt.subset2,
                                       Route==this.route)
            this.pred.subset3 <- subset(this.pred.subset2,
                                        Route==this.route)
            for (this.dose in unique(this.cvt.subset3$Dose))
              if (this.dose %in% tolower(this.pred.subset3$Dose))
            {
              this.cvt.subset4 <- subset(this.cvt.subset3,
                                         Dose==this.dose)              
              this.pred.subset4 <- subset(this.pred.subset3,
                                          Dose==this.dose)
              obs.times <- signif(sort(unique(this.cvt.subset4$Time)),4)
              this.cvt.subset4means <- NULL
              for (this.time in obs.times)
              {
                this.cvt.subset5 <- subset(this.cvt.subset4,
                                           signif(Time,4)==this.time)
                new.row <- data.frame(
                  Compound=this.compound,
                  DTXSID=this.dtxsid,
                  CAS=this.cas,
                  Species=this.species,
                  Route=this.route,
                  Dose=this.dose,
                  Time=this.time,
                  Conc.pred=this.pred.subset4[obs.times == this.time,"Conc_est"],
                  stringsAsFactors=FALSE)
                new.row <- merge(new.row,
                                 this.cvt.subset5[,c(
                                   "CAS","Media","Value","Source","calc_loq")], by="CAS")
                colnames(new.row)[colnames(new.row)=="Value"] <- "Conc.obs"
                cvt.table <- rbind(cvt.table,new.row)
                this.cvt.subset4means <- rbind(this.cvt.subset4means,
                                           data.frame(
                                             Time=this.time,
                                             Value=mean(this.cvt.subset5$Value,
                                                        na.rm=TRUE)))
              }
              Cmax.obs <- max(this.cvt.subset4means$Value, na.rm=TRUE)
              Cmax.pred <- max(this.pred.subset4[,
                                                 "Conc_est"], 
                               na.rm=TRUE)
              if (length(obs.times)>1)
              {
                AUC.obs <- AUC(this.cvt.subset4means$Time,
                               this.cvt.subset4means$Value)
                if (length(obs.times) != 
                    length(this.pred.subset4[
                      !duplicated(this.pred.subset4[,"Time"]),
                      "Conc_est"])) browser()
                AUC.pred <- AUC(obs.times,
                                this.pred.subset4[
                                  !duplicated(this.pred.subset4[,"Time"]),
                                                  "Conc_est"])
              } else {
                AUC.obs <- NA
                AUC.pred <- NA
              }
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
                stringsAsFactors=FALSE)
              stats.table <- rbind(stats.table,new.row)            
            }
          }
        
      }
    
  }
  
  cvt.table$QSPR <- label
  stats.table$QSPR <- label
  return(list(
    cvt=cvt.table,
    stats=stats.table))
}


# This function does the level III statistic comparisons:
maketkstatpreds <- function(
  CvT.data,
  cvtfits,
  label)
{
  chems.good.1comp <- get_cheminfo(model="1compartment") 
  out.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
    if (this.cas %in% cvtfits$CAS)
  {
    if (this.cas %in% chems.good.1comp)
    {
      this.subset1 <- subset(CvT.data,CAS==this.cas)
      this.compound <- this.subset1$Compound[1]
      this.dtxsid <- this.subset1$DTXSID[1]
      this.subset1 <- subset(cvtfits,CAS==this.cas)
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
  out.table$QSPR <- label
  return(out.table)
}
