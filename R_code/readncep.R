require(matlab)
ncep=read.table("ncep_temp.txt",skip=5,header = TRUE)
ncep=ncep[,c(1,3)]
ncep$sinlat=sin(ncep$lat/180*pi)
ncep$hemsyn=(ncep$Temperature + flipud(ncep$Temperature))/2
