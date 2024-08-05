# This function does the level II concentration comparisons:
makeCvTpreds <- function(CvT.data,label,model.args)
{
  nonvol.chems <- suppressWarnings(get_cheminfo(model="pbtk"))
  vol.chems <- suppressWarnings(get_cheminfo(model="gas_pbtk"))
  
  cvt.table <- NULL
  stats.table <- NULL
  for (this.cas in unique(CvT.data$CAS))
  if (this.cas %in% c(nonvol.chems, vol.chems))
  {
    this.model <- NULL
    #if (this.cas %in% nonvol.chems) this.model <- "solve_pbtk"
    #else if (this.cas %in% vol.chems) this.model <- "solve_gas_pbtk"
    # Just use gas_pbtk for everything:
    this.model <- "solve_gas_pbtk"
    
    if (!is.null(this.model))
    {
      print(paste(label,": ",this.cas,sep=""))
      this.subset1 <- subset(CvT.data,CAS==this.cas & 
        !is.na(Time_Days) & 
        !is.na(Conc_mgpL))
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
            obs.times <- signif(sort(unique(this.subset4$Time_Days)),4)

            if (this.model=="solve_pbtk")
            {
              params <- suppressWarnings(do.call("parameterize_pbtk",
                                                 args=c(list(
                                                   ,
                                                   species=this.species,
                                                   default.to.human=TRUE),
                                                   model.args)))
            } else if (this.model=="solve_gas_pbtk")
            {
              params <- suppressWarnings(do.call("parameterize_gas_pbtk",
                                                 args=c(list(
                                                   chem.cas=this.cas,
                                                   species=this.species,
                                                   default.to.human=TRUE),
                                                   model.args)))
            }
            if ("Caco2.options" %in% names(model.args))
            {
              if ("keepit100" %in% names(model.args[["Caco2.options"]]))
              {
                if (model.args[["Caco2.options"]][["keepit100"]])
                  params$Fabsgut <- 1
              }
            }
            pred <- suppressWarnings(do.call(this.model,
              args=c(
                list(parameters=params,
                     times=sort(unique(c(seq(0,2,0.05),obs.times))),
                     species=this.species,
                     iv.dose=(this.route=="iv"),
                     dose=this.dose,
                     default.to.human=TRUE,
                     suppress.messages=TRUE,
                     input.units='mg/kg',
                     output.units = "mg/L",
                     exp.conc=0),
                model.args)))
# Units of pred are Cplasma: mg/L and time: days
            pred <- subset(pred,pred[,"time"]>0.0001)
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
              this.subset5 <- subset(this.subset4,signif(Time_Days,4)==this.time)
  #            print(this.subset5)
              new.row <- data.frame(
                Compound=this.compound,
                DTXSID=this.dtxsid,
                CAS=this.cas,
                Species=this.species,
                Route=this.route,
                Dose=this.dose,
                Time_Days=this.time,
                Conc.pred=pred[pred[,"time"] == this.time,"Cven"],
                stringsAsFactors=F)
              new.row <- merge(new.row,
                this.subset5[,c(
                  "CAS","Media","Conc_mgpL","Source","calc_loq")], by="CAS")
              colnames(new.row)[colnames(new.row)=="Conc_mgpL"] <- "Conc.obs"
              if (any(tolower(new.row[,"Media"])=="blood"))
              {
                new.row[tolower(new.row[,"Media"])=="blood","Conc.pred"] <-
                  signif(
                    Rb2p * 
                      new.row[tolower(new.row[,"Media"])=="blood","Conc.pred"],
                    4)
              }
              cvt.table <- rbind(cvt.table,new.row)
              this.subset4means <- rbind(this.subset4means,
                data.frame(
                  Time_Days=this.time,
                  Conc_mgpL=mean(this.subset5$Conc_mgpL,na.rm=0)))
            }
            Cmax.obs <- max(this.subset4means$Conc_mgpL,na.rm=T)
            Cmax.pred <- max(pred[,"Cven"],na.rm=T)
            if (length(unique(this.subset4means$Time_Days))>1)
            {
              AUC.obs <- AUC(this.subset4means$Time_Days,this.subset4means$Conc_mgpL)
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
  
  out <- calc_cvt_stats(cvt.table, stats.table)
  
  return(out)
}

calc_cvt_stats <- function(cvt.table, stats.table)
{
  # Set a minimal value for predictions and observations using limit of quantification (LOQ)
  # We treat all "low" values as the same, where we define low with loq:
  for (this.col in c("Conc.pred","Conc.obs"))
  {
    cvt.table[cvt.table[,this.col] < cvt.table$calc_loq, this.col] <- 
      cvt.table[cvt.table[,this.col] < cvt.table$calc_loq, "calc_loq"]
  }
  min.loq <- min(cvt.table[,"calc_loq"],na.rm=TRUE)
  for (this.col in c("Cmax.obs","Cmax.pred","AUC.obs","AUC.pred"))
  {
      which.rows <- stats.table[,this.col] < min.loq
      which.rows[is.na(which.rows)] <- FALSE
      stats.table[which.rows, this.col] <- min.loq
  }
  for (this.row in 1:dim(stats.table)[1])
  {
    stats.table[this.row, "Cmax.RMSLE"] <- calc_RMSLE(stats.table[this.row, ],
                                                      obs.col="Cmax.obs",
                                                      pred.col="Cmax.pred")
    stats.table[this.row, "AUC.RMSLE"] <- calc_RMSLE(stats.table[this.row, ],
                                                     obs.col="AUC.obs",
                                                     pred.col="AUC.pred")
  }
  
  # Calculate absolute fold error for all observations
  cvt.table$AbsFE <- calc_AbsFE(cvt.table)
  
  # Calculate chemical-specific RMSLE and AAFE:
  rmsle <- aafe <- rmsle.early <- rmsle.late <- aafe.early <- aafe.late <- list()
  for (this.chem in unique(cvt.table$DTXSID))
  {
    this.data <- subset(cvt.table,DTXSID==this.chem)
    rmsle[[this.chem]] <- calc_RMSLE(this.data, zeroval = -Inf)
    aafe[[this.chem]] <- calc_AAFE(this.data)
    
    # Separate data into early and late times:
    this.data.early <- this.data.late <- NULL
    for (this.study in unique(this.data$Source))
    {
      this.sourcedata <- subset(this.data,Source==this.study)
      mid.time <- median(unique(this.sourcedata$Time_Days))
      this.data.early <- rbind(this.data.early,
                               subset(this.sourcedata,
                                      Time_Days < mid.time))
      this.data.late <- rbind(this.data.late,
                              subset(this.sourcedata,
                                     Time_Days >= mid.time))
    }
    rmsle.early[[this.chem]] <- calc_RMSLE(this.data.early, zeroval = -Inf)
    aafe.early[[this.chem]] <- calc_AAFE(this.data.early)
    rmsle.late[[this.chem]] <- calc_RMSLE(this.data.late, zeroval = -Inf)
    aafe.late[[this.chem]] <- calc_AAFE(this.data.late)
  }
  
  return(list(
    cvt=cvt.table,
    stats=stats.table,
    rmsle=rmsle,
    aafe=aafe,
    rmsle.early=rmsle.early,
    rmsle.late=rmsle.late,
    aafe.early = aafe.early,
    aafe.late = aafe.late
    ))
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
                             !is.na(Time_Days) & 
                             !is.na(Conc_mgpL))
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
              obs.times <- signif(sort(unique(this.subset4$Time_Days)),4)
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
                this.subset5 <- subset(this.subset4,signif(Time_Days,4)==this.time)
                new.row <- data.frame(
                  Compound=this.compound,
                  DTXSID=this.dtxsid,
                  CAS=this.cas,
                  Species=this.species,
                  Route=this.route,
                  Dose=this.dose,
                  Time_Days=this.time,
                  Conc.pred=pred[obs.times == this.time],
                  stringsAsFactors=FALSE)
                new.row <- merge(new.row,
                                 this.subset5[,c(
                                   "CAS",
                                   "Media",
                                   "Conc_mgpL",
                                   "Source",
                                   "calc_loq")], by="CAS")
                colnames(new.row)[colnames(new.row)=="Conc_mgpL"] <- "Conc.obs"
                cvt.table <- rbind(cvt.table,new.row)
                this.subset4means <- rbind(this.subset4means,
                                           data.frame(
                                             Time_Days=this.time,
                                             Conc_mgpL=mean(this.subset5$Conc_mgpL,na.rm=0)))
              }
              Cmax.obs <- max(this.subset4means$Conc_mgpL,na.rm=T)
              Cmax.pred <- max(pred,na.rm=T)
            #  if (Cmax.pred == 0) browser()
              if (length(obs.times)>1)
              {
                AUC.obs <- signif(AUC(this.subset4means$Time_Days,
                                      this.subset4means$Conc_mgpL), 4)
                AUC.pred <- signif(AUC(obs.times, pred), 4)
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
  
  out <- calc_cvt_stats(cvt.table, stats.table)
  
  return(out)
}

# This function does the level III statistic comparisons:
maketkstatpreds <- function(
  CvT.data,
  cvtfits,
  label)
{
  chems.good.1comp <- suppressWarnings(get_cheminfo(model="1compartment"))
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
        cl.pred <- 1/suppressWarnings(calc_total_clearance(
          chem.cas=this.cas,
          species=this.species,
          default.to.human=TRUE,
          suppress.messages=TRUE))
        ke.pred <- cl.pred/vd.pred 
        thalf.pred <- signif(log(2)/ke.pred,3)
        new.tab <- data.frame(
          Compound=this.compound,
          DTXSID=this.dtxsid,
          CAS=this.cas,
          Species=this.species,
          Vd.obs=vd.obs,
          Vd.pred = vd.pred,
          thalf.obs=thalf.obs,
          thalf.pred = thalf.pred,
          cl.pred = cl.pred,
          stringsAsFactors=F)
  #      new.tab <- merge(new.tab,this.subset2[,
  #        c("CAS","Reference")],by="CAS")
        out.table <- rbind(out.table,new.tab)
      }
    }
  }
  out.table$QSPR <- label
  return(out.table)
}

## calculate absolute fold error
calc_AbsFE <- function(level2tab)
{
  #Absolute fold error (AbsFE) : abs(log10(pred/obs))
  AbsFE<-signif(abs(log10(level2tab$Conc.pred/level2tab$Conc.obs)),
                          4)

  #Turn -Inf and Inf into finite value
  AbsFE[AbsFE=="Inf"]<-4   ######## lots because lots of $Conc.pred is zero

  #Turn NaN into 0 , because the NaN is a result of matching 0/0...and hence it was correctly identified.
  AbsFE[AbsFE=="NaN"]<-0  #none

  return(AbsFE)
}


## Absolute Average Fold Error (AAFE) :  
# 10^((1/n)*sum(abs(FE)))
# use this. It's good to see fold error
calc_AAFE <- function(in.table,
                      AbsFE.col="AbsFE",
                      sigfig=4)
{
    return(signif(
    10^(mean(in.table[,AbsFE.col],
                  na.rm=TRUE)),
    sigfig))
}

## Root mean squared log10 Error (RMSLE): 
# sqrt(mean(log10(Xpred+1)-log10(Xobs+1))2)
# use this
calc_RMSLE <- function(in.table,
                       pred.col="Conc.pred",
                       obs.col="Conc.obs",
                       sigfig=4,zeroval=1e-6)
{
# Get red of NA's:
  in.table <- subset(in.table, 
                      !is.na(in.table[,pred.col]) & 
                      !is.na(in.table[,obs.col]))

# How to handle very small values which diverge to negative infinity:
  for (this.col in c(pred.col, obs.col))
    in.table[in.table[,this.col] < zeroval,this.col] <- zeroval

  return(signif(
    sqrt(mean((log10(in.table[,pred.col] /
                       in.table[,obs.col]))^2,
                   na.rm=TRUE)),
    sigfig))
}

## Median relative predictive error (MPRE)
calc_MRPE <- function(level2.table)
{
  return(signif(
    median(level2.table$RPE,
                na.rm=TRUE),
    4))
}

## Relative predictive error (MPRE)
# (pred - obs) / obs
calc_RPE <- function(obs, pred)
{
  return((pred - obs) / obs)
}


# Function for formatting tick labels:
scientific_10 <- function(x) {                                  
  out <- gsub("1e", "10^", scientific_format()(x))              
  out <- gsub("\\+","",out)                                     
  out <- gsub("10\\^01","10",out)                               
  out <- parse(text=gsub("10\\^00","1",out))                    
}  

# Generate standard table comparing predictions and observations for multiple
# chemical lists
makestatstable <- function(obspred.table,
                           obs.col,
                           pred.col,
                           chem.lists, # Different lists for subsetting the chemicals
                           logrsquared = TRUE,
                           average.per.chemical = FALSE)
{
  stats.table <- data.frame()
  for (this.qspr in unique(obspred.table$QSPR))
    for (list.index in 1:length(chem.lists))
    {
      this.chem.list <- chem.lists[[list.index]]
      this.chem.list.name <- names(chem.lists)[list.index]
      these.data <- subset(obspred.table, 
                           QSPR == this.qspr & 
                             DTXSID %in% this.chem.list)
      these.data <- subset(these.data,
                           !is.na(these.data[,pred.col]) &
                           !is.na( these.data[,pred.col]))
      if (!average.per.chemical)
      {
        this.count <- length(unique(these.data$DTXSID))
        this.nobs <- dim(these.data)[1]
        if (this.nobs > 1)
        {
          if (logrsquared)
          {
            this.fit <- lm(log10(eval(these.data[,obs.col])) ~ eval(log10(these.data[,pred.col])))
          } else {
            this.fit <- lm(eval(these.data[,obs.col]) ~ eval(these.data[,pred.col]))
          }
          this.rsquared <- signif(summary(this.fit)$r.squared, 2)
        } else {
          this.rsquared <- NA
        }
        this.rmsle <- signif(calc_RMSLE(these.data, 
                                        pred.col = pred.col,
                                        obs.col = obs.col), 2)
        stats.table[paste("RMSLE",this.chem.list.name), this.qspr] <- paste(
          this.rmsle, " (",
          this.count, ", ",
          this.nobs, ")",
          sep="")
        stats.table[paste("RSquared",this.chem.list.name), this.qspr] <- paste(
          this.rsquared, " (",
          this.count, ", ",
          this.nobs, ")",
          sep="")
      } else {
        mean.rmsle <- NULL
        for (this.chem.index in 1:length(unique(these.data$DTXSID)))
        {
          this.chem <- unique(these.data$DTXSID)[this.chem.index]
          mean.rmsle[this.chem.index] <- signif(calc_RMSLE(subset(these.data,
                                                                  DTXSID == this.chem), 
                                                           pred.col = pred.col,
                                                           obs.col = obs.col), 3)
        }
        this.mean.rmsle <- signif(mean(mean.rmsle),2)
        this.count <- length(mean.rmsle)
        stats.table[paste("RMSLE",this.chem.list.name), this.qspr] <- paste(
          this.mean.rmsle, " (",
          this.count, ")",
          sep="")
      } 
    } 
  return(stats.table)
}

makestatstable2 <- function(this.table,
                            stats.list,
                            chem.lists,
                            label="",
                            sigfig=3)
{
  for (this.qspr in names(stats.list))
    for (this.list in names(chem.lists))
    {
      these.stats <- stats.list[[this.qspr]][
        names(stats.list[[this.qspr]]) %in% chem.lists[[this.list]]]
      this.count <- length(these.stats)
      this.mean <- signif(mean(unlist(these.stats), na.rm=TRUE), sigfig)
      this.table[paste(this.list, label), this.qspr] <- paste(
        this.mean, " (",
        this.count, ")",
        sep="")
    }
  return(this.table)
}


kelim_substitution <- function(data,
                               invivo,
                               Vd_data,
                               clint_column,
                               fup_column) {
  
  stopifnot(is.data.frame(data),
            is.data.frame(invivo),
            is.data.frame(Vd_data),
            is.character(clint_column),
            is.character(fup_column))
  
  if (!all(c("CAS", clint_column, fup_column) %in% names(data))) {
    stop("Error: CAS or specified clint/fup columns not in data.frame")
  }
  if (!all(c("CAS", "kelim") %in% names(invivo))) {
    stop("Error: CAS or kelim columns not in fit table")
  }
  if (!all(c("CAS", "Vd.pred") %in% names(Vd_data))) {
    stop("Error: CAS or Vd.pred columns not in Vd_data")
  }
  
  # Inner join by CAS number to add kelim column for ease of calculation
  data <- merge(data, invivo[c("CAS", "kelim")])
  data <- merge(data, Vd_data[c("CAS", "Vd.pred")])
  data$`invivoAdjusted.Fup` <- NA
  data$`invivoAdjusted.Clint` <- NA
  
  for (y in seq_len(nrow(data))) {
    
    this_data <- data[y, ]
    this_cas <- this_data[["CAS"]]
    # Message to see which chemical the process is on.
    print(this_cas)
    
    this_paramset <- suppressWarnings(httk::parameterize_pbtk(chem.cas = this_cas))
    Qgfr <- this_paramset[["Qgfrc"]]
    Qli <- this_paramset[["Qliverf"]]
    Rb2p <- this_paramset[["Rblood2plasma"]]
    Vdist <- this_data["Vd.pred"]
    kelim <- this_data["kelim"]
    
    # Summarizing some variables (left hand terms)
    lht <- kelim * Vdist - Qgfr
    
    # here we calculate and create two independent columns for 
    # adjusted CLint and fup
    fup_data <- this_data[fup_column]
    clint = 1/(fup_data * (lht - (1/(Qli * Rb2p))))
    clint <- signif(clint, 4)
    
    clint_data <- this_data[fup_column]
    fup = 1/(clint_data * (lht - (1/(Qli * Rb2p))))
    fup <- signif(fup, 4)
    
    # Before writing, there needs to be checks for negative values
    # I will set either to NA if below zero (or if fup > 1 set to 1)
    clint <- ifelse(clint < 0, NA, clint)
    fup <- ifelse(fup < 0, NA, fup)
    fup <- ifelse(fup > 1, 1, fup)
    
    
    data[y, "invivoAdjusted.Fup"] <- fup
    data[y, "invivoAdjusted.Clint"] <- clint
    # When using these values, always use *_column for the other value in
    # subsequent analyses
    
  }
  
  
  return(data)
  
  
}
