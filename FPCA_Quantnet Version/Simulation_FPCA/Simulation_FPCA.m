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
            [Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simulate(N,j*N*T+1,T) ;            
            
            %%estimate the variances  
            sigma=zeros(N,1);
            parfor(i=1:N)
                sigma(i)=variance(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
            end
            
            %%estimation using three different methods see describtion of fpca1_gpu and individual
            [hX2b, V2b, loadsb,Meansmo2bb,Db]   =fpca_gpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma,0);
            [hX2a, V2a, loadsa,Meansmo2ba,Da]   =fpca_gpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,2,'cpu',sigma,0);
            [hX2i]                              =individual(Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma);

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

