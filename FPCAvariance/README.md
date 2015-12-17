
[<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/banner.png" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **FPCAvariance** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet: FPCAvariance

Published in: 'Functional Principal Component Analysis for Derivatives of
High-Dimensional Spatial Curves.'

Description: 'Calculates two-dimensional variance estimator, based on Hall and Marron (1990).'

Keywords: simulation, empirical, local polynomial, surface, nonparametric, derivative

See also: 'FPCAgpu, FPCASimulation, FPCAmultiloc, FPCAsimulate_input, FPCAindividual, FPCAepan'

Author: Heiko Wagner

Submitted:  Maria Grith


Input:
 - X1 : x1 coordinate
 - X2 : x2 coordinate
 - Y  : function value


Output:
 - sigma : variance estimator

```


```matlab
function [sigma] = FPCAvariance( X1, X2, Y)

x1     =X1;
x2     =X2;
kernel ='';
N1     = length( X1 );
T      = length( x1 );
d      =2;
h0     =T^(-2/(4+d));
Wges   =zeros(N1,T);

for (j=1:T)
    Xmx1 = X1 - x1(j);
    Xmx2 = X2 - x2(j);
    if(strcmp(kernel,'Gauss')==1)
        W = normpdf(Xmx1/h0) .*normpdf(Xmx2/h0)  ;
        %'using gaussian kernel'
    else
        W = FPCAepan(Xmx1/h0) .*FPCAepan(Xmx2/h0)  ;
    end
    Wges(:,j) =W/sum(W);
end

v     = T - 2 * sum(diag(Wges)) + sum(sum( Wges.^2)) ;
sigma = 1/v* sum( (Y-Wges'*Y).^2 );
end

```
