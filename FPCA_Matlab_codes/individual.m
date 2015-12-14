% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: Low dimensional decomposition of derivatives for surfaces
%              via method M(d) with automatic bandwith selection
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
%L          - Number of Dimensions
%Fc         - Triscattered intepolated curves
%x1minc     - Min x1 value
%x2minc     - Min x2 value
%x1maxc     - Max x1 value
%x2maxc     - Max x2 value
%cgridx     - Joined Grid x1
%cgridy     - Joined Grid x2
%c5unil     - Moneyness Axis
%c2unil     - Maturity Axis
%c3unil     - Obervationw with error
%c3unilr    - Obervations wo. error
%N          - Number of Days         
%Tmon       - Obervations  Monetary axis
%Tmat       - Obervations Matuiry axis
  
% ------------------------------------------------------------------------------ 
% Output:      
%Low dimensional decomposition using L Dimensions
%Eigenvectors of second derivative
%Corresponding loadings
%Mean curve
%Eigenvalues

% ------------------------------------------------------------------------------ 
% Keywords:    FPCA, Surfaces, Derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2015/08/13 
% ------------------------------------------------------------------------------ 

function [hX2r,mx,my] = individual(Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,method,comp,sigma)
mxe=mx;
mye=my;
%%grid construction begin           
%%to compute M we need a random grid
x1min   =max( x1minc );
x2min   =max( x2minc );
x1max   =min( x1maxc );
x2max   =min( x2maxc );
Tint=length(cell2mat(c5unil(1)) );
my      =x1min + ( x1max-x1min ).*rand( Tint , 1 );     
mx      =x2min + ( x2max-x2min ).*rand( Tint , 1 );
densy   =x2max-x2min;
densx   =x1max-x1min;

%%fit the data to a common
Yint    =[];
Tmat    =zeros(N,1);
Tmon    =zeros(N,1);

scores=[];
parfor i=1:N
    
    xachs   =cell2mat(  c2unil(i) );
    yachs   =cell2mat(  c5unil(i) );
    Lp      =[xachs.^0, xachs, xachs.^2, xachs.^3, xachs.^4, xachs.^5, yachs.^1, yachs.^2, yachs.^3, yachs.^4];
    Vxre    =cell2mat(  c3unil(i) );
    scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
    Tmat(i) =length(unique(cell2mat(  c2unil(i) ) ));    
    Tmon(i) =length(unique(cell2mat(  c5unil(i) ) ));
    ddum    =Fc{i}( mx(:) , my(:) );
    Yint    =[Yint ddum(:)];
end
X               =Yint- repmat( mean( Yint , 2 ), [1 N]);
X( isnan(X) )   =0;

%%grid construction end
%%ADD derivative constant for 1st derivative HERE

Cdp=1.0006;   %%constant for derivatives taken from Fan 1996  

%%estimate sigma
if(sigma==0 ) 
    sigma=zeros(N,1);
    
    parfor i=1:N
        sigma(i)=varaince(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
    end
end
sigmae  =sigma;

%%%estimate h0 in each direction
h1a =[];
h2a =[];

Lppd4x=[mx*0, mx*0, mx*0, 0*6*mx.^0, 24*mx.^0, 120*mx.^1, my*0, my*0, my*0, 24*my*0];
Lppd4y=[mx*0, mx*0, mx*0, 6*mx*0, 24*mx*0, 120*mx*0, my*0, my*0, 0*6*my.^0, 24*my.^0];
regsx=Lppd4x*scores; 
regsy=Lppd4y*scores;

parfor i=1:N 
    h1a= [h1a (1/Tmat(i))^(1/10)* (sigmae(i)/( densx* (mean(regsy(:,i) ).^2 ) )  )^(1/10)]; 
    h2a= [h2a (1/Tmon(i))^(1/10)* (sigmae(i)/( densy* (mean(regsx(:,i) ).^2    ) )  )^(1/10)];
end

Xsmoa   =[] ;
Xsmo2a  =[] ;
parfor i=1:N 
    if comp=='gpu'
        [XmiS0a XmiS1a XmiS2a XmiS3a]=multiloc( gpuArray(cell2mat(  c5unil(i) )),gpuArray(cell2mat(  c2unil(i) )) , gpuArray(cell2mat(  c3unil(i) )),mye ,mxe ,h1a(i), Cdp *h2a(i) ,3,'Gauss' ); %%estimate loc poly -----------switched to lower case
    else
        [XmiS0a XmiS1a XmiS2a XmiS3a]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) ,cell2mat(  c3unil(i) ),mye ,mxe ,h1a(i), Cdp *h2a(i) ,3,'Gauss' ); %%estimate loc poly -----------switched to lower case
    end
    Xsmoa   =[Xsmoa XmiS0a'];
    Xsmo2a  =[Xsmo2a XmiS2a'];
end
if comp=='gpa'
    Xsmoa=gpuArray(Xsmoa);
    Xsmo2a=gpuArray(Xsmo2a);  
end

%%final estimator   
hX2r        =Xsmo2a;
end

%plotting options  
%scatter3(my,mx,hX2r(:,2))
%scatter3(my,mx,realD(:,2))
%scatter3(my,mx, Xsmo2a(:,2))

