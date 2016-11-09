% clear variables and close windows
clear all;
close all;
clc;

%%add your path here; genpath(folderName) returns a path string that includes folderName and multiple levels of subfolders below folderName.
%addpath(genpath('your path'));

%%% Start READ REAL DATA

% read expiration dates
fid = fopen('ExpirationDates.txt','rt');
indata = textscan(fid, '%s %s', 'HeaderLines',1);
fclose(fid);
expdata = indata{1,1};
[ye, me, de] = datevec(expdata,'dd-mmm-yy');

% read DAX index
[a1,a2]=textread('dax_index.txt','%s %f');
dax=horzcat(a2);
date_dax=datenum(a1,'dd.mm.yyyy');

% read VDAX index
[a1v,a2v]=textread('vdax.txt','%s %f');
vdax=horzcat(a2v);
date_vdax=datenum(a1v,'dd.mm.yyyy');

% read interest rate
[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16]=textread('IRs.dat','%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
date_ir=datenum(b1,'dd.mm.yyyy');
IR=horzcat(b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

% read call data
[c11,c12,c13,c14,c15,c16,c17]=textread('C_2002.txt', '%s %s %s  %f %f %f %f');
[c21,c22,c23,c24,c25,c26,c27]=textread('C_2003.txt', '%s %s %s  %f %f %f %f');
[c31,c32,c33,c34,c35,c36,c37]=textread('C_2004.txt', '%s %s %s  %f %f %f %f');
[c41,c42,c43,c44,c45,c46,c47]=textread('C_2005.txt', '%s %s %s  %f %f %f %f');
[c51,c52,c53,c54,c55,c56,c57]=textread('C_2006.txt', '%s %s %s  %f %f %f %f');
[c61,c62,c63,c64,c65,c66,c67]=textread('C_2007.txt', '%s %s %s  %f %f %f %f');
[c71,c72,c73,c74,c75,c76,c77]=textread('C_2008.txt', '%s %s %s  %f %f %f %f');
[c81,c82,c83,c84,c85,c86,c87]=textread('C_2009.txt', '%s %s %s  %f %f %f %f');
[c91,c92,c93,c94,c95,c96,c97]=textread('C_2010.txt', '%s %s %s  %f %f %f %f');
[c101,c102,c103,c104,c105,c106,c107]=textread('C_2011.txt', '%s %s %s  %f %f %f %f');

c1g=[c11;c21;c31;c41;c51;c61;c71;c81;c91;c101];
c2g=[c12;c22;c32;c42;c52;c62;c72;c82;c92;c102];
c3g=[c13;c23;c33;c43;c53;c63;c73;c83;c93;c103];
c4g=[c14;c24;c34;c44;c54;c64;c74;c84;c94;c104];
c5g=[c15;c25;c35;c45;c55;c65;c75;c85;c95;c105];
c6g=[c16;c26;c36;c46;c56;c66;c76;c86;c96;c106];
c7g=[c17;c27;c37;c47;c57;c67;c77;c87;c97;c107];

% define maturies for the interest rates
T=[7/365 14/365 21/365 1/12 2/12 3/12 4/12 5/12 6/12 7/12 8/12 9/12 10/12 11/12 1];

days=unique(c2g);
day=datenum(c2g, 'dd/mm/yyyy');
N=length(days);

%initiate container to store the results
vdaxc=[];
odax=[];
R=[];
daysa=[];

cgridx=[];
cgridy=[];

x1minc=[];
x2minc=[];
x1maxc=[];
x2maxc=[];

for d=2:(N+1)
    tday=datenum(days((d-1)), 'dd/mm/yyyy');
    select= find( day==tday ) ;

    c1=c1g(select);
    c2=c2g(select);
    c3=c3g(select);
    c4=c4g(select);
    c5=c5g(select);
    c6=c6g(select);
    c7=c7g(select);

    refdmy_c=datenum(c2, 'dd/mm/yyyy');
    cd=[refdmy_c, c4, c5, c6, c7];
    numdate=unique(refdmy_c);

    % compute maturities
    clear tauc
    clear aa

    for i=1:length(cd)
        aa(i)=find(ye==c5(i) & me==c4(i));
        tauc(i)=(datenum(expdata(aa(i)),'dd-mmm-yy')-refdmy_c(i,1))/365;
    end

    % find stock price
    dday=find(date_dax==numdate);
    s=dax(dday);
    odax=[odax,s];

    if(s>0)
        daysa=[daysa, dday ];

        % find volatility index
        vdday=find(date_vdax==numdate);
        sv=vdax(vdday);
        vdaxc=[vdaxc,sv];

        % find interest rate
        rday=find(date_ir==numdate);
        R=IR(rday,:);

        % interpolate interest rate
        rc=interp1(T,R,tauc,'linear','extrap');

        % date rescaled_strike call_price interest_rate maturity
        cdd=[cd(:,1) cd(:,4)./(s*exp(rc'/100.*tauc')) cd(:,5)./(s*exp(rc'/100.*tauc')) rc'/100 tauc'];

        % remove observations with maturity larger than one year
        cdd(find(tauc>1),:)=[];

        % sort maturities
        [cdd(:,5), Im] = sort( cdd(:,5) );
        cdd(:,2)=cdd(Im,2);
        cdd(:,3)=cdd(Im,3);

        % sorting moneyness
        matu=unique(cdd(:,5));
        for l=1:length(matu)
            Ik=find(cdd(:,5)==matu(l) );
            [cdd(Ik,2),Il] =sort(cdd(Ik,2));
            cdd(Ik,3)=cdd((min(Ik)-1)+Il,3);
        end

        % interpolate data for FPCA
        F = TriScatteredInterp(cdd(:,5),cdd(:,2),cdd(:,3), 'linear');
        % store interpolated curves
        Fc{d-1}=F;


        % find smallest common grid
        x1minc=[ x1minc min(cdd(:,5))];
        x2minc=[ x2minc  min(cdd(:,2))];
        x1maxc=[ x1maxc max(cdd(:,5))];
        x2maxc=[ x2maxc max(cdd(:,2))];

        cgridx=union(cgridx,cdd(:,5));
        cgridy=union(cgridy,cdd(:,2));
        c5unil{d-1}= cdd(:,5);   % maturity
        c2unil{d-1}=cdd(:,2);    % moneyness
        c3unil{d-1}=cdd(:,3);    % call price
    end
end
%%% End READ REAL DATA

%%% Start ESTIMATION


x1min   =median( x1minc );
x2min   =median( x2minc );
x1max   =median( x1maxc );
x2max   =1.4; % or use median( x2maxc );

% construct common grid
my      =x1min + ( x1max-x1min ).*rand( 512 , 1 ); % maturity
mx      =x2min + ( x2max-x2min ).*rand( 512 , 1 ); % moneyness

% estimate the variances
sigma=zeros(N,1);
parfor i=1:N
    sigma(i)=FPCAvariance(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
end

L=3; % low dimensional space parameter
method=1; % decomposition of the dual covariance matrix of the observed curves
% estimation using FPCAgpu
[hX2b, V2b, loadsb,Meansmo2bb,Db]=FPCAgpu(L,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unil,N,mx,my,method,'cpu',sigma,1);

%%% End ESTIMATION
