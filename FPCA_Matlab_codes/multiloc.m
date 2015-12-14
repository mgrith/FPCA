% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: 2 Dimensional local polynomial estimator 
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
%X1 - x1 coordinate, 
%X2 - x2 coordinate,  
%Y  - function value
%x1 - x1 output coordinate
%x2 - x2 output coordinate 
%h1 - bandwith in x1 direction
%h2 - bandwith in x2 direction
%p  - degree of polynomial
 
% ------------------------------------------------------------------------------ 
% Output:      
%fit0  - Smoothed curves
%fit1  - Smoothed 1. derivative
%fit2  - Smoothed 2. derivative
%fit3  - Smoothed 3. derivative
%h1    - bandwith in x1 direction
%h2    - bandwith in x2 direction

% ------------------------------------------------------------------------------ 
% Keywords:    local polynomial surface estimator, derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2015/12/08 
% ------------------------------------------------------------------------------ 


function [fit0, fit1, fit2, fit3, W0, W1, W2, W3] = multiloc( X1, X2,  Y, x1,x2, h1,h2,p,kernel)

N1 = length( X1 );
Na = length( x1 );

function [fit0, fit1, fit2, fit3, W0, W1, W2, W3] = getpoint(j)
Xmx1 = X1 - x1(j);
Xmx2 = X2 - x2(j);   
Xn   = [Xmx1(:) Xmx2(:)]; 
Z    = ones( N1, 2*p+1 ); 
for m=1:p
   Z(:,(2*m):(2*m+1) ) = Xn.^m;
end

if(strcmp(kernel,'Gauss')==1) 
    W=diag( normpdf(Xmx1/h1)/h1 .*normpdf(Xmx2/h2)/h2  );
    %'Using Gaussian Kernel'
else
    W=diag( epan(Xmx1/h1)/h1 .*epan(Xmx2/h2)/h2  );
end

A    =Z' * W * Z;    
Ws   = pinv(A) * Z' * W;
b    = Ws * Y;
fit1 =0;
W1   = 0;
fit2 = 0;
W2   = 0;
fit3 = 0;
W3   = 0;
  
%deriv=0;
fit0 =   b(1,1);
W0   =   sum(Ws(1,:).^2) ;
  if(p>1)     
      %deriv=1;
      fit1 =factorial(1)*b(3,1);
      W1   = factorial(1)^2*sum(Ws(3,:).^2);
  end
  if(p>2)    
      % deriv=2;
      fit2 =factorial(2) * b(5,1);  
      W2   = factorial(2)^2*sum(Ws(5,:).^2);
      %deriv=3;
      fit3 =factorial(3)*b(7,1);  
      W3   = factorial(2)^2*sum(Ws(7,:).^2);
  end
end
  
[fit0, fit1, fit2, fit3, W0, W1, W2, W3]=arrayfun(@getpoint,1:Na); 
end

