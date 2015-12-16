
[<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/banner.png" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **FPCAepan** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet: FPCAepan

Published in: Functional Principal Component Analysis for Derivatives of High-Dimensional Spatial Curves

Description: Epanechnikov kernel function

Keywords: Epanechnikov kernel, smoothing, nonparametric

See also: FPCAgpu, FPCAsimulation, FPCAmultiloc, FPCAsimulate_input, FPCAindividual, FPCAvariance 

Author: Heiko Wagner

Submitted:  Maria Grith


Input:   
- t : time points


Output:   
- k : kernel function

```


```matlab
function [k]=FPCAepan(t);
k= 3/4*(1-t.^2);
k((abs(t)-1)>0)=0;
end

```
