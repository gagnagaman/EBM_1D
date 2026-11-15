library(matlab)
source("ebm.R")
source("readncep.R")

# default values
out=ebm()

with(out,plot(lat,Temperature))
with(ncep,lines(lat,hemsyn))

with(out,plot(sinlat,Temperature))
with(ncep,lines(sinlat,hemsyn))

# play with D
out=ebm(D=0.77)
with(out,plot(lat,Temperature,type="l"))
with(ncep,points(lat,hemsyn))

# play with A
out=ebm(A=180)
with(out,plot(lat,Temperature,type="l"))
with(ncep,points(lat,hemsyn))

# play with B
out=ebm(B=1.6)
with(out,plot(lat,Temperature,type="l"))
with(ncep,points(lat,hemsyn))

