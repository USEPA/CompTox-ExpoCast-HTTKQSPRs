# Function to clear all TK data out of httk:
clear_httk <- function(
  target.dtxsids=NULL, 
  target.env=.GlobalEnv,
  species=c("Human","Rat"))
{
  for (this.species in species)
  {
    delete.table <- chem.physical_and_invitro.data[,c("Compound","CAS","DTXSID")]

    if (!is.null(target.dtxsids)) delete.table <- subset(delete.table,
      DTXSID %in% target.dtxsids)
    
    if (length(delete.table[,"DTXSID"])>0)
    {
      these.indices <- chem.physical_and_invitro.data[,"DTXSID"] %in%
        delete.table[,"DTXSID"]
      
      for (this.col in c(
        "Clint","Funbound.plasma","Rblood2plasma","Clint.pValue"))
      {
        col.name <- paste(this.species,this.col,sep=".")
        chem.physical_and_invitro.data[these.indices, col.name] <- NA
        chem.physical_and_invitro.data[these.indices, paste(
          col.name,"Reference",sep=".")] <- "Deleted"
      }
    }
  }
    
  assign("chem.physical_and_invitro.data", 
    chem.physical_and_invitro.data, envir=target.env)
}                                  
