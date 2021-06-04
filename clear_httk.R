# Function to clear all TK data out of httk:
clear_httk <- function(
  target.dtxsids=NULL, 
  target.env=.GlobalEnv,
  species=c("Human","Rat"))
{
  for (this.species in species)
  {
    delete.table <- get_cheminfo(info="all",species=this.species)
    if (!is.null(target.dtxsids)) delete.table <- subset(delete.table,
      DTXSID %in% target.dtxsids)
    if (dim(delete.table)[1]>0)
    {
      delete.table$Clint <- NA
      delete.table$Funbound.plasma <- NA
      delete.table$Rblood2plasma <- NA
      delete.table$Clint.pValue <- NA
      chem.physical_and_invitro.data <- 
        add_chemtable(delete.table,
          current.table=chem.physical_and_invitro.data,
          data.list=list(
          Compound="Compound",
          CAS="CAS",
          DTXSID="DTXSID",
          Funbound.plasma="Funbound.plasma",
          Rblood2plasma="Rblood2plasma",
          Clint.pValue="Clint.pValue",
          Clint="Clint"),
          species=this.species,
          reference="Deleted",
          overwrite=T,
          allow.na=T)      
    }
  }
    
  assign("chem.physical_and_invitro.data", 
    chem.physical_and_invitro.data, envir=target.env)
}                                  
