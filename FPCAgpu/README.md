
[<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/banner.png" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **FPCAgpu** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet: FPCAgpu

Published in: Functional Principal Component Analysis for Derivatives of High-Dimensional Spatial Curves

Description: 'Computes low dimensional decomposition of derivatives for surfaces using automatic bandwidth selection.'

Keywords: simulation, empirical, FPCA, surface, derivative, density

See also: FPCASimulation, FPCAepan, FPCAmultiloc, FPCAsimulate_input, FPCAindividual, FPCAvariance

Author: Heiko Wagner

Submitted:  Maria Grith

Input:
            - L          : Number of Dimensions
            - Fc         : Triscattered interpolated curves
            - x1minc     : Min x1 value
            - x2minc     : Min x2 value
            - x1maxc     : Max x1 value
            - x2maxc     : Max x2 value
            - cgridx     : Joined Grid x1
            - cgridy     : Joined Grid x2
            - c5unil     : Moneyness Axis
            - c2unil     : Maturity Axis
            - c3unil     : Observations with error
            - c3unilr    : Observations wo. error
            - N          : Number of Days
            - mx         : Output Grid Monetary axis
            - my      	 : Output Grid Maturity axis
            - method     : Computation method (1=M0, 2=Md)
            - comp       : GPU or No GPU (experimental, GPU usually much slower)
            - sigma      : Estimated or ture sigma
            - app        : Application (1=Yes, 0=No) effects diagonal correction

Output:
            - hX2r       : Low dimensional decomposition using L Dimensions
            - V2a        : Eigenvectors of second derivative
            - loadsa     : Corresponding loadings
            - Meansmo2b  : Mean curve
            - Da         : Eigenvalues
            - mx         : Output Grid Monetary axis
            - my         : Output Grid Maturity axis

```


```matlab
function [hX2r, V2a, loadsa,Meansmo2b,Da,mx,my] = FPCAgpu(L,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,method,comp,sigma,app)
L   =L-1;
mxe =mx;
mye =my;

%%grid construction begin
%%to compute M we need a random grid
x1min   =max( x1minc );
x2min   =max( x2minc );
x1max   =min( x1maxc );
x2max   =min( x2maxc );

Tint    =max(length(cell2mat(c5unil(1)) ),512);
my      =x1min + ( x1max-x1min ).*rand( Tint , 1 );
mx      =x2min + ( x2max-x2min ).*rand( Tint , 1 );
densy   =x2max-x2min;
densx   =x1max-x1min;

%Fit the data to a common grid
Tmat    =zeros(N,1);
Tmon    =zeros(N,1);
Yint    =[];
parfor i=1:N
    Tmat(i) =length(unique(cell2mat(  c2unil(i) ) ));
    Tmon(i) =length(unique(cell2mat(  c5unil(i) ) ));
    ddum    =Fc{i}( mx(:) , my(:) );
    Yint    =[Yint ddum(:)];
end
X               =Yint- repmat( mean( Yint , 2 ), [1 N]);
X( isnan(X) )   =0;

%%grid construction end
Cdp=1.0006;   %%constant for derivatives taken from Fan 1996

%%estimate hat(sigma)
if(sigma==0 )
    sigma=zeros(N,1);
    parfor i=1:N
        sigma(i)=varaince(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
    end
end
sigmae  =sigma;
Lppd4x  =[mx*0, mx*0, mx*0, 0*6*mx.^0, 24*mx.^0, 120*mx.^1, my*0, my*0, my*0, 24*my*0];
Lppd4y  =[mx*0, mx*0, mx*0, 6*mx*0, 24*mx*0, 120*mx*0, my*0, my*0, 0*6*my.^0, 24*my.^0];
h1m     =[];
h2m     =[];
scores  =[];
parfor i=1:N

    %%estimate optimal individual smoothing parameter
    %%Construct polynomial estimator for h1
    xachs   =cell2mat(  c2unil(i) );
    yachs   =cell2mat(  c5unil(i) );
    Lp      =[xachs.^0, xachs, xachs.^2, xachs.^3, xachs.^4, xachs.^5, yachs.^1, yachs.^2, yachs.^3, yachs.^4];
    Vxre    =cell2mat(  c3unil(i) );
    scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];

    h2i=(max(cell2mat(  c2unil(i) ))-min(cell2mat(  c2unil(i) )))   *Tmat(i);
    h1i=(max(cell2mat(  c5unil(i) ))-min(cell2mat(  c5unil(i) )))   *Tmon(i);
    h1m     =[h1m h1i];
    h2m     =[h2m h2i];
end

W       =[];
Xsmo2m  =[];

for i=1:N

    if(method==1)
        if comp=='gpu'
            [XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=FPCAmultiloc( gpuArray(cell2mat(  c5unil(i) )),gpuArray(cell2mat(  c2unil(i) )) , gpuArray(cell2mat(  c3unil(i) )),my ,mx ,h1m(i)^(-1/3),h2m(i)^(-1/3),1,'Epan'  ); %%estimate loc poly
        else
            [XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=FPCAmultiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,h1m(i)^(-1/3),h2m(i)^(-1/3),1,'Epan'  ); %%estimate loc poly
        end

        Xsmo2m =[Xsmo2m XmiS0m'] ;
        W      = [W, W0];
    else
        if comp=='gpu'
            [XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=FPCAmultiloc( gpuArray(cell2mat(  c5unil(i) )),gpuArray(cell2mat(  c2unil(i) )) , gpuArray(cell2mat(  c3unil(i) )),my ,mx ,h1m(i)^(-1/13),h2m(i)^(-1/13),7,'Epan'  ); %%estimate loc poly
        else
            [XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=FPCAmultiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,h1m(i)^(-1/10),h2m(i)^(-1/10),7,'Epan'  ); %%estimate loc poly
        end
        Xsmo2m =[Xsmo2m XmiS2m'] ;
        W      = [W, W2];
    end

end

if comp=='gpa'
    Xsmo2m=gpuArray(Xsmo2m);
end

Xa            =Xsmo2m;
Muc           =(1/size(Xa,1))*(Xa'*Xa - diag(mean(W).*sigmae)) ;
M             =Muc - (1/size(Xa,1))*Xa'*repmat(mean(Xsmo2m,2),[1 N]) - (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*Xa + (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*repmat(mean(Xsmo2m,2),[1 N]);
[Vduala, Da]  =eig(M);
[Da,Ind]      =sort(Da(Da>0),'descend');
Da2           =Da
Vduala        =Vduala(:,Da>0);
Vduala        =Vduala(:,Ind);
Da            =diag(Da);
L             = min( (size(Da2,1)-1),L );
regsx         =Lppd4x*scores;
regsy         =Lppd4y*scores;
h1a           =zeros(1,(L+1));
h2a           =zeros(1,(L+1));

for i=1:(L+1)
    h1a(i) = (1/mean(Tmat) )^(1/10)* (mean(sigmae)/( densx*mean( (regsy*Vduala(:,i)).^2  )))^(1/10);
    h2a(i) = (1/mean(Tmon) )^(1/10)* (mean(sigmae)/( densy*mean( (regsx*Vduala(:,i)).^2  )))^(1/10);
end
h1a     =ones(1,N)*mean(h1a);
h2a     =ones(1,N)*mean(h2a);
Xsmoa   =[] ;
Xsmo2a  =[] ;
parfor i=1:N
    if comp=='gpu'
        [XmiS0a XmiS1a XmiS2a XmiS3a]=FPCAmultiloc( gpuArray(cell2mat(  c5unil(i) )),gpuArray(cell2mat(  c2unil(i) )) , gpuArray(cell2mat(  c3unil(i) )),mye ,mxe ,h1a(i), Cdp *h2a(i) ,3,'Gauss' ); %%estimate loc poly
    else
        [XmiS0a XmiS1a XmiS2a XmiS3a]=FPCAmultiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) ,cell2mat(  c3unil(i) ),mye ,mxe ,h1a(i), Cdp *h2a(i) ,3,'Gauss' ); %%estimate loc poly
    end
    Xsmoa   =[Xsmoa XmiS0a'];
    Xsmo2a  =[Xsmo2a XmiS2a'];
end
if comp=='gpa'
    Xsmoa  =gpuArray(Xsmoa);
    Xsmo2a =gpuArray(Xsmo2a);
end

Meansmo2b   =mean(Xsmo2a')';                %mean curve
V2a         =(Xsmo2a )*Vduala*sqrt(Da)^-1;	%eigenvectors 2nd derivative smoothed estimate
loadsa      = sqrt(Da)*Vduala';             %loadings

%%final estimator
hX2r        =repmat(Meansmo2b,[1 N])+ V2a(:,1:L+1)*loadsa(1:L+1,:);
end

```
