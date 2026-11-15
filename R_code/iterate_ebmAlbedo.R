iterate_ebmAlbedo <- function(Temp,So,disp=TRUE) {
  # Using a temperature profile Temp that was
  # generated with a call to ebm or ebmAlbedo
  # and a value for the solar constant, iterate
  # the ebmAlbedo model to figure out what is the
  # equilibrium solution for both temperature and
  # albedo given the So.
   
  # By default it prints the status of work
  # This can be turned off.
  
  diffoff=1
  Talb=Temp
  counter=0
  while (diffoff >1e-5) {
   newout=ebmAlbedo(Talb,So) 
   Talb_new=newout$Temperature
   diffoff=sum(abs(Talb_new-Talb))
   Talb=Talb_new
   counter = counter +1
   if (counter > 50) break
   if (disp) print(paste("solar constant",So," Counter",counter,
                         " difference ",diffoff,
                         "average albedo",mean(newout$albedo)))
  }
  return(newout)
}
