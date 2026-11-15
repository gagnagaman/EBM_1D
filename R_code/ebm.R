ebm <- function(N=NULL,A=NULL,B=NULL,Diff=NULL) {
  # Usage:
  #  T=ebm() 
  # or 
  #  T=ebm(N=200,A=204,B=1.9,Diff=0.42)
  # or some other values.
  #
  # A Code that sets up a steady state version of the zonally
  # averaged energy balance model of North (1975, JAS 32) and 
  # solves for T(x) where x=sin(lat).
  #
  # The model equation is:
  #
  # Qsw+Q_T-Qlw =0
  #
  # where Qsw is the solar insolation
  # Qlw is the heat radiation to space 
  # and Q_T is the net heat flux divergence.
  #
  # Approximate formulas for Qsw are
  # 
  # Qsw=So/4* 1/4+(5-3x^2)*(1-albedo(x)) (North)
  # A similar one was used by Fanning and Weaver (see below)
  # Qsw=So/4*(1 -0.492*(3x^2-1)/2)*(1-albedo(x)) 
  # which can be written as 
  # Qsw=1370/4*(1 -0.492*P2(x))*(1-albedo(x)) 
  # where So is the solar constant (1370 W/m^2) and P2(x) is the
  # second legendre polynomial. The albedo is often specified as a
  # scalar, roughly 0.3 - 0.4 on a global but can  also be specified,
  # either as a function of  temperature, or as a function of
  # latitude. 
  #
  # An example of a latitude function can be found in Fanning and
  # Weaver 1996, JGR vol 101)  albedo(x)=0.315+0.3615*x^2.
  # This formulation has an average of 0.43,  with equatorial
  # values of 0.315 rising to 0.68 at the poles.
  # Temperature dependent albedo values are often used for the
  # exploration of the small ice cap instability. An  example is to
  # assume ice cover (and albedo of 0.6) if T<-10 C but to use 0.33
  # otherwise. This code must be modified slightly to allow for
  # temperature dependend albedo.  
  #
  # The formula for Qlw is the standard linearized version
  # Qlw=A+B*T where A,B are constants and T is temperature in deg. C.
  # A range of  values can be found in the litterature for of A &
  # B. Budyko's original values were A=211 and B=1.9, but values
  # often seen are   A=203.3 W/m^2 and B=2.09 W/m^2C.
  # See Graves et.al (1993, JGR vol 98) for details.
  #
  # Q_T is the divergence of the heatflux. In North's model the heat
  # flux is parameterized as down gradient diffusion
  # F=D*dT/dx. Using spherical coordinates to calculate divergences
  # and gradients, and then using x=sin(lat) one obtains
  # 
  # Q_T= d/dx [ D(x) (1 - x^2) dT/dx ]
  #
  # where D(x) is the effective diffusivity constant, and the 1-x^2
  # part takes care of the fact that there is less area to average
  # over closer to the poles.
  #
  # D can either be a scalar or a function of latitude. Scalar values
  # for D are  found by tuning D so that calculated T fits
  # observations. Values ranging from  0.35 to 0.65, depending on the
  # value of other coefficients. The lower range being preferred by North. 
  #
  # D can also be a function of latitude. For instance Lindzen and
  # Farrel (JAS, 1977) suggest a value of
  #              D(x)=0.45*[1+9*exp(-(x/sin(30))^6)]
  # to account for the fact that the Hadley cell is more efficient at 
  # heat transfer than mid-latitude systems.
  
  # Using constant coefficients for D & albedo this code has been
  # compared against an analytical solution (in terms of legendre
  # polynomials, P0 and P2) and the maximum difference was found to
  # be less than  1.8E-5 degC.
  #
  # Halldor Bjornsson 04/2005
  # Translated from matlab in August 2018. 
  # Sorry for the lack of R-ish code.
  
  #####################################################
  #       Load library & define helper function
  #####################################################
  library(matlab)
  legendre2 <- function(x) (3*x^2 -1)/2
  #####################################################
  #         If not defined in function call...
  #####################################################
  if (is.null(N)) N=200;          # number of gridcells: 
                                  #   how good is your computer?
  if (is.null(Diff)) Diff=0.42;   # for diffusion, code can 
                                  #   handle latitudinally varying diff.
  if (is.null(A)) A=204;          # Sources vary on number...
  if (is.null(B)) B=1.9;          # Sources vary on number...
  
  #####################################################
  #              Define scalar parameters
  #####################################################
  So=1370;     # Solar constant 
  alpha=0.33;  # for albedo. code can handle a 
               #   latitudinally varying alpha
  
  #########################################################
  # make grid. We solve for one hemisphere (x btw 0,1) but 
  # could just as well use both hemispheres (x btw -1 and 1).
  #########################################################
  
  xe=linspace(0,1,N);      # cell-edges, this is where the fluxes are defined
  xm=xe[1:(N-1)]+diff(xe)/2; # mid cells, where Qsw and Qlw  (and T)
                           # are defined.
  dx=xe[2]-xe[1];
  
  
  ##############################################################
  #              Define latitudinally varying parameters
  ##############################################################
  
  # Setup insolation albedo and solar heating
  p=legendre2(xm);
  S=So/4*(1-0.492*p);
  albedo=ones(size(p))*alpha; # modify this if latitudinally
  # varying albedo is preferred.
  # have tried: albedo=(0.315+0.3615*xm.^2)/1.2 
  # in R the albedo vector must be wrapped in a matrix:
  # albedo=matrix((0.315+0.3615*xm^2)/1.2,nrow=1)
  
  Qsw=(1-albedo)*S;
  
  # Set up diffusion
  D=Diff*ones(1,(N-2)); # modify this if latitudinally varying
  # diffusion is preferred.
  # An example could be Lindzen and Farrel
  # D=matrix(0.45*(1+9*exp(-1*(xe[2:(length(xe)-1)]/sin(30/180*pi))^6)),
  #         nrow=1) 
  ##############################################################
  #              Set up Matrices
  ##############################################################
  
  # Set up a matrix for B*T
  Bmat=eye(N-1)*B;
  
  # make the matrix wich takes care of the "(1-x^2)*D(x)" bit
  diagline=D * (1-xe[2:(N-1)]^2);
  Dmat=zeros(N-2,N-2);
  for (i in c(1:(N-2)))  Dmat[i,i]=diagline[i];
  
  # set matrix that calculates dT/dx at xe(2:N-1). Note that this is
  # done in terms of grid cell averges, hence this is an exact
  # formulation (not O(1)). 
  
  FT=zeros(N-2,N-1);
  for (i in c(1:(N-2))) {
    FT[i,i]=-1;
    FT[i,i+1]=1;
  }
  
  FT=FT/dx;
  
  # set up matrix to calculate d/dx (Flux) 
  # based on cell averaged differening. O(x^3) accuracy. Our life is
  # made easier by the fact that  F(1)=F(N) = 0. 
  QF=zeros(N-1,N-2);
  for (i in c(1:(N-2))) {
    QF[i,i]=1;
    QF[i+1,i]=-1;
  }
  QF=QF/dx;
  
  ##########################################
  # Run the model....
  #########################################
  # We have d/dx((1-x^2)*D* d/dx T) - B*T. 
  # The first d/dx is QF, the second is FT
  # The (1-x^2)*D is in Dmat, and B is in Bmat.
  # Hence in matrix form this is
  # QF*(Dmat*FT*T)-Bmat*T or (QF*Dmat*FT - Bmat)*T
  
  # set up matrix 
  
  M=QF%*%Dmat%*%FT - Bmat;
  
  # solve the model (anti-climatic isn't it?)
  Temperature=solve(M)%*%t(A-Qsw);
  
  #variables to return
  OLR=A+B*Temperature
  
  outvar=list(sinlat=xm,
              lat=asin(xm)*180/pi,
              Temperature=as.numeric(Temperature),
              Qsw=as.numeric(Qsw),
              Qlw=as.numeric(OLR) )
 
  return(outvar)
}
