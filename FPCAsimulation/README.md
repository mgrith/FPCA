
[<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/banner.png" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **FPCAsimulation** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet: FPCAsimulation

Published in: 'Functional Principal Component Analysis for Derivatives of
High-Dimensional Spatial Curves'

Description: 'Starts a simulation for varying number of curves and observations per curve for repeated samples,
and performs a comparison of the methods.'

Keywords: simulation, FPCA, surface, derivative, density

See also: FPCAgpu, FPCAepan, FPCAmultiloc, FPCAsimulate_input, FPCAindividual, FPCAvariance

Author: Heiko Wagner

Submitted:  Maria Grith


Input:
 - reps : Repetitions +1
 - Nn   : Vector of to simulate curves per trail
 - Tt   : Vector of to simulate observations per curve

Output:
 - meanM : Matrix of Mean MSE with entries for each value of Nn and Tt
 - varV  : Matrix of variances of MSE with entries for each value of Nn and Tt
 - medM  : Matrix of medians of MSE with entries for each value of Nn and Tt
 - iqrM  : Matrix of IQR of MSE with entries for each value of Nn and Tt

```


```matlab
% clear variables and close windows
clear all;
close all;
clc;

%%add your path here; genpath(folderName) returns a path string that includes folderName and multiple levels of subfolders below folderName.
%addpath(genpath('your path'));
%%start simulation
%initiate container to store the results
meanM =[];
varV  =[];
medM  =[];
iqrM  =[];

%setting fixed variables
reps  =2;
Nn    =[10 25];
Tt    =[50 250];

for (N=Nn)
    for (T=Tt)

        hXm1b =[];
        hXm1a =[];
        hXm1i =[];
        j     =1;

        while(j<reps)

            %%simulate data
            [Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=FPCAsimulate_input(N,j*N*T+1,T) ;

            %%estimate the variances
            sigma=zeros(N,1);
            parfor(i=1:N)
                sigma(i)=FPCAvariance(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
            end

            %%estimation using three different methods see describtion of FPCAgpu and individual
            [hX2b, V2b, loadsb,Meansmo2bb,Db]   =FPCAgpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma,0);
            [hX2a, V2a, loadsa,Meansmo2ba,Da]   =FPCAgpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,2,'cpu',sigma,0);
            [hX2i]                              =FPCAindividual(Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma);

            %%store MSE
            hXm1b =[hXm1b,  mean( (hX2b-realD).^2) ];
            hXm1a =[hXm1a,  mean( (hX2a-realD).^2) ];
            hXm1i =[hXm1i,  mean( (hX2i-realD).^2) ];
            j     =j+1
        end

        mean(hXm1b)   %M_0 Method
        mean(hXm1a)   %M_d Method
        mean(hXm1i)   %Individual

        var(hXm1b)   %M_0 Method
        var(hXm1a)   %M_d Method
        var(hXm1i)   %Individual

        median(hXm1b)   %M_0 Method
        median(hXm1a)   %M_d Method
        median(hXm1i)   %Individual

        iqr(hXm1b)   %M_0 Method
        iqr(hXm1a)   %M_d Method
        iqr(hXm1i)   %Individual

        meanM =[meanM , [mean(hXm1b) mean(hXm1a) mean(hXm1i)]'  ]
        varV  = [varV , [var(hXm1b) var(hXm1a) var(hXm1i)]' ]
        medM  = [medM , [median(hXm1b) median(hXm1a) median(hXm1i)]' ]
        iqrM  = [iqrM , [iqr(hXm1b) iqr(hXm1a) iqr(hXm1i)]' ]
    end
end

```
