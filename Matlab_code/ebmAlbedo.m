 function [T,xm]=ebmAlbedo(T,So)
% calculate the temperature that results from adjusting the albedo
% acording to the temperature. Take in So so it can be varied for 
% studies of the small ice cap instability.

% Halldor Bjornsson 04/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Define scalar parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(T)+1;
%N=200;      % number of gridcells: how good is your computer?
%So=1370 ;   % Solar constant 
Diff=0.42;   % for diffusion, code can handle latitudinally varying diff.
alpha=0.33; % for albedo. code can handle a latitudinally varying alpha
A=204;      % Sources vary on number...
B=1.9;      % Sources vary on number...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make grid. We solve for one hemisphere (x btw 0,1) but 
% could just as well use both hemispheres (x btw -1 and 1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xe=linspace(0,1,N);      % edges, this is where the fluxes are defined
xm=xe(1:N-1)+diff(xe)/2; % mid cells, where Qsw and Qlw  (and T)
                         % are defined.
dx=xe(2)-xe(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Define latitudinally varying parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup insolation albedo and solar heating
p=legendre(2,xm);
S=So/4*(1-0.492*p(1,:));
albedo=ones(size(p(1,:)))*alpha; % modify this if latitudinally
                                 % varying albedo is preferred.
				 % have tried:
                                 % albedo=(0.315+0.3615*xm.^2)/1.2;
albedo(T<-10)=0.6;               % if temperature is below -10 we
                                 % assume ice cover with mean
                                 % albedo of 0.6!			 
				
Qsw=(1-albedo)'.*S';

% Set up diffusion
D=Diff*ones(1,N-2); % modify this if latitudinally varying
		    % diffusion is preferred.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Set up Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up a matrix for B*T
Bmat=eye(N-1)*B;
 
% make the matrix wich takes care of the "(1-x^2)*D(x)" bit
diagline=D .* (1-xe(2:N-1).^2);
Dmat=zeros(N-2,N-2);
for i=1:N-2
  Dmat(i,i)=diagline(i);
end
		    
% set matrix that calculates dT/dx at xe(2:N-1). Note that this is
% done in terms of grid cell averges, hence this is an exact
% formulation (not O(1)). 

FT=zeros(N-2,N-1);
for i=1:N-2
  FT(i,i)=-1;
  FT(i,i+1)=1;
end
FT=FT/dx;

% set up matrix to calculate d/dx (Flux) 
% based on cell averaged differening. O(x^3) accuracy. Our life is
% made easier by the fact that  F(1)=F(N) = 0. 
QF=zeros(N-1,N-2);
for i=1:N-2
QF(i,i)=1;
QF(i+1,i)=-1;
end
QF=QF/dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the model....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We have d/dx((1-x^2)*D* d/dx T) - B*T. 
% The first d/dx is QF, the second is FT
% The (1-x^2)*D is in Dmat, and B is in Bmat.
% Hence in matrix form this is
% QF*(Dmat*FT*T)-Bmat*T or (QF*Dmat*FT - Bmat)*T

% set up matrix 

M=QF*Dmat*FT - Bmat;

% solve the model (anti-climatic isn't it?)

T=inv(M)*(A-Qsw);

return  %(uncomment return if function statement at top is
% uncommented). 


