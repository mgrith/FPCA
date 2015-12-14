function [Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simulate(N,k,T)   
  
x1minc  =[];
x2minc  =[];
x1maxc  =[];
x2maxc  =[];
cgridx  =[];  
cgridy  =[];
Qs      =[];
Xs      =[];
taugrid =[];

%simulate data
randn('state',k+6);
ino  =randn(1,N+2)';

r     =0.02;  %interest rate 
sigma =0.1;   %observation error
Tout  =256;   %number of observations true model

%simulate stock price
st(1)=100;  
for i=2:(N+1)
    st(i)=100;
end

%parameters of the mixtures 
m   =[0.4 0.7 0.1]; %mean 
s   =[0.5 0.3 0.3]; %standard deviation

%weights of the mixture
rand('state',k+1);
p =abs( randn(3,N+1 )  );
p =p./repmat(sum(p),3,1);

for i=2:(N+1)
    Xj   =[];
    Yj   =[];
    Cj   =[];
    qj   =[];
    tauj =[];
    
    for j=1:T
        tau = 0.2+rand( 1 , 1 )*0.5 ; %maturity
        x   = 0.5+rand( 1 , 1 )*1.3;  %moneyness
        
        for l=1:3
            qu(l,:) =1./(x*sqrt(2*pi*s(l)^2*tau)).*exp(-1/2*((log(x)- (m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau))).^2);  %density components of the mixture
            Ca(l,:) =exp(m(l)*tau)*normcdf((log(1/x)+(m(l)+s(l)^2/2)*tau)/(s(l)*sqrt(tau)))-x.*normcdf((log(1/x)+(m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau))); %components of the call price
        end
        
        Ca(isnan(Ca)) = 0 ;
        X(j,:,i)      =x;                              %observation Grid K
        Q(j,:,i)      =p(:,i)'*qu;                     %value of risk neutral density
        C(j,:,i)      =p(:,i)'*Ca;                     %call price
        Y(j,:,i)      =C(j,:,i)+sigma*randn(1,1 );     %obervation with error
        Xj            =[Xj ; X(j,:,i)'];               %strike
        Yj            =[Yj ; Y(j,:,i)'];               %calls
        Cj            =[Cj ; C(j,:,i)'];               %calls w/o error
        tauj          =[tauj ; tau*ones(size(C,2),1)]; %moneyness
        qj            =[qj ; Q(j,:,i)'];
    end
    
    cdd=[i*ones(length(Xj),1) Xj Yj r*ones(length(Xj),1) tauj ]; % Day strike call_price interest_date maturity
    
    %%simulated completed
    %store interpolated curves
    F           =TriScatteredInterp(cdd(:,2),cdd(:,5),cdd(:,3), 'nearest');
    Fc{i-1}     =F;
    
    %find smallest common grid and scale prices
    x1minc       =[x1minc min(cdd(:,5))];        %maturity
    x2minc       =[x2minc  min(cdd(:,2))];       %moneyness
    x1maxc       =[x1maxc max(cdd(:,5))];
    x2maxc       =[x2maxc max(cdd(:,2))];
    cgridx       =union(cgridx,cdd(:,5));
    cgridy       =union(cgridy,cdd(:,2));
    c5unil{i-1}  =cdd(:,5);                      %maturity
    c2unil{i-1}  =cdd(:,2);                      %strike/s*exp(rtau)
    c3unil{i-1}  =cdd(:,3);                      %call/stock_price
    c3unilr{i-1} =Cj/st(i).*exp(r*tauj);
end

%store true curves
x1min =min( x1minc );
x2min =min( x2minc );
x1max =max( x1maxc );
x2max =max( x2maxc );
my    =0.2+rand( Tout , 1 )*0.5;
mx    =0.5+rand( Tout , 1 )*1.3;        
X     =mx;
T     =my;

for l=1:3
    QQ(l,:) =1./(X.*sqrt(2*pi*s(l)^2.*T)).*exp(-1/2*((log(X)- (m(l)-s(l)^2/2)*T)./(s(l)*sqrt(T))).^2);    % equivalent of qu for X and T
end

realD=QQ'*p(:,2:end);

end
