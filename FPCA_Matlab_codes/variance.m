% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: 2 Dimensional Variance estimator, Based on Hall, Marron 1990
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
%X1 - x1 coordinate, 
%X2 - x2 coordinate,  
%Y  - function value
 
% ------------------------------------------------------------------------------ 
% Output:      
%sigma

% ------------------------------------------------------------------------------ 
% Keywords:    local polynomial surface estimator, derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2015/12/08 
% ------------------------------------------------------------------------------ 


function [sigma] = variance( X1, X2,  Y)
x1 =X1;
x2 =X2;

kernel='';

N1   = length( X1 );
T    = length( x1 );
d    =2;
h0   =T^(-2/(4+d));
Wges =zeros(N1,T);

for (j=1:T)  
    Xmx1 = X1 - x1(j);
    Xmx2 = X2 - x2(j);
    if(strcmp(kernel,'Gauss')==1)
        W = normpdf(Xmx1/h0) .*normpdf(Xmx2/h0)  ;
        %'using gaussian kernel'
    else
        W = epan(Xmx1/h0) .*epan(Xmx2/h0)  ;
    end
    Wges(:,j) =W/sum(W);
end
  
v     = T - 2 * sum(diag(Wges)) + sum(sum( Wges.^2)) ; 
sigma = 1/v* sum( (Y-Wges'*Y).^2 );
end


