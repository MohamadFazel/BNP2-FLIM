function K=calCOV(X1,X2,T,L)
%calCOV calculate the squared exponential correlation matrix
%
%INPUT:
%   X1: First set of coordinates
%   x2: Second set of coordinates
%   T:  Parameter of the correlation function
%   L:  Parameter of the correlation function
%
%OUTPUT:
%   K:  Correlation matrix between the first and second sets of coordinates
%
%Created by:
%   Mohamadreza Fazel, (Prsse lab, 2022)   
%

Dist = pdist2(X1,X2);
K =(T^2).*exp(-Dist.^2/L^2/2);
end