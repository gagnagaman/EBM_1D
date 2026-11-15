analyT = function(slat=seq(0,1,by=0.001),A=NULL,B=NULL,Diff=NULL) {
  #   The analytical function derived for the ebm
  if (is.null(Diff)) Diff=0.42;   # Sources vary on number...
  if (is.null(A)) A=204;          # Sources vary on number...
  if (is.null(B)) B=1.9;          # Sources vary on number...
  legendre2 <- function(x) (3*x^2 -1)/2
  
  Q=342 
  a=0.7
  p=legendre2(slat);
  beta = -0.492
  Qsw=Q*a*(1 + beta * p);
  Temperature = (Q*a - A)/B  + beta * Q*a/(B + 6 *Diff) * p
  
 outvar=list(
   slat=slat,
   Temperature=Temperature,
   Qsw=Qsw,
   Qlw=A+B*Temperature)
 
 return(outvar)
}

 

 