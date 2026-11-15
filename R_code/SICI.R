# Run the Small Ice Cap Instability experiment
# code translated from matlab.

library(matlab)
source("ebmAlbedo.R")
source("ebm.R")
source("iterate_ebmAlbedo.R")

So_start=1370  # Initial value of the solar constant. The ebm
               # will not ice over with this value.

# Lets go upwards in solar constant and examine the ice free branch
SoUpp=linspace(So_start,So_start*1.30,40) # So values to use for
                                          # examining the ice free branch
Tupp=zeros(size(SoUpp)); # vector to store the temperatures
					  
outAlb=ebm() # Make an initial temperature profile to use


print('Examine the upper part of the ice free branch')

for (i in (1:length(SoUpp))) { 
  outAlb=iterate_ebmAlbedo(outAlb$Temperature,SoUpp[i],disp=FALSE)
  Tupp[i]=mean(outAlb$Temperature)
}

# Examine the lower part of the ice-free branch and and lower So until
# the earth becomes ice covered. I have already determined that So=1277
# is enough to trigger a snowball state. So if I go down in So by 15%
# that will be more than enough. 
SoDown=linspace(So_start,So_start*(1-0.15),71)
Tdown=zeros(size(SoDown))

print('Examine the rest of the ice free branch down to the ice covered one')

outAlb=ebm() # Again an initial temperature profile to use

for (i in (1:length(SoDown))) {
  outAlb=iterate_ebmAlbedo(outAlb$Temperature,SoDown[i],disp=FALSE)
  Tdown[i]=mean(outAlb$Temperature)
}

# Examine the ice covered branch from So=1250 with increasing
# So until the earth goes back to an ice free state. 

So_cold=1150;
outAlb=iterate_ebmAlbedo(outAlb$Temperature,So_cold,disp=FALSE)

SoUpp2=linspace(So_cold,So_start*1.30,250);
Tupp2=zeros(size(SoUpp2));

print('examine the ice covered state increasing So until ice free state')

for (i in (1:length(SoUpp2))) {
  outAlb=iterate_ebmAlbedo(outAlb$Temperature,SoUpp2[i],disp=FALSE)
  Tupp2[i]=mean(outAlb$Temperature)
}

# plot results
plot(SoUpp2,Tupp2,type="l",
     xlab="Solar constant (W/m2)",ylab="Global temperature (°C)")
lines(SoUpp,Tupp,col=2)
lines(SoDown,Tdown,col=3)

# plot 2
plot(SoUpp2,Tupp2,cex=0.1,pch=19,
     xlab="Solar constant (W/m2)",ylab="Global temperature (°C)")
points(SoUpp,Tupp,col=2,cex=0.1)
points(SoDown,Tdown,col=3,cex=0.1)
