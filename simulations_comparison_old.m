% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: Starts a simulation as described in the article and perform 
%              both in the article described methodes
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
         %Number of days
         %Risk-free interest rate
         %Drift
         %Volatility
         %Observation error
         %Number of Dimensions 
         %Obervations Monetary axis
         %Obervations Matuiry axis
% Additional Inputs Estimation: 
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
% Keywords:    FPCA, Surface, 2nd dervative
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2015/05/28 
% ------------------------------------------------------------------------------ 

clear all 
%addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_neu')); 
addpath(genpath('I:\OwnCloud2\Uni\Programs\Quantanet Version_neu')); 
addpath(genpath('I:\OwnCloud2\Uni\Programs\Quantanet Version_GPU\')); 

addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_GPU'));  


addpath(genpath('C:\Users\Heiko Wagner\Desktop\Quantanet Version_GPU')); 


%%START SIMULATION
%%INPUTS: 
         %Number of days
         %We observe calls & stock price every month
         %Risk-free interest rate
         %Drift
         %Volatility
         %observation error
         %Number of Dimensions +1  ! Just to confuse people 
         %Obervations Monetary axis
         %Obervations Matuiry axis
 
%N=10;
%T=50;

meanM=[]
varV=[]
medM=[]
iqrM=[]



m=0
for (N=[10 25])
for (T=[50 250])

m=m+1    

hXm1b=[];
hXm1a=[];
hXm1i=[];
j=1;

N=15
T=100
[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simu1(N, 1/12, 0.02, 0.03, 0.18, 0.1, j*N*T+1,T, 256 ) ;    %One dimensional        
%[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simu2(N, 1/12, 0.02, 0.03, 0.18, 0.1, j*N*T+1,T, 256 ) ;    

cell2mat( c5unil(1) ) %<--- grid fixed at 0.5 ! can be dropped 
scatter( cell2mat(c2unil(1) ) ,cell2mat(c3unil(1) ) )
scatter3( cell2mat(c5unil(1) ) ,cell2mat(c2unil(1) ) ,cell2mat(c3unil(1) ) )
% 
% %%Perform M�ller algorithm
% 
% addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_GPU\PACE_matlab-master\release2.17\PACE')); 
% 
% addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_GPU\PACE_matlab-master\release2.17\PACE\PACE-DER')); 
% 
% p = setDerOptions('nder',0:2)      
% [yy] = FPCder(c3unil, c2unil ,p);
% ypred = getVal(yy,'y_pred');
% y0 = ypred{1};     %estimated curve
% y1 = ypred{2};     %estimated first derivative
% y2 = ypred{3};     %estimated second derivative
% 
% 
% scatter( cell2mat(c2unil(1) ) ,cell2mat(c3unil(1) ) )
% hold on
% plot(cell2mat( yy(20) ),cell2mat(y0(1)) )
% 
% 
% plot(cell2mat( yy(20) ),cell2mat(y2(1)) )

%%Estimate the variances once to save time
sigma=zeros(N,1);
parfor(i=1:N)
sigma(i)=varaince( [  cell2mat(  c2unil(i) ) ],cell2mat(  c3unil(i) ) );
end

Xc=[]
for i=1:N
    Xc{i}= [ cell2mat(c5unil(i)),cell2mat(c2unil(i)) ]
end


[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(2,Fc,5,[min(x1minc),min(x2minc)],[max(x1maxc),max(x2maxc)],Xc,c3unil,c3unilr,[mx,my],1,sigma);
scatter3(x(:,1),x(:,2),hX2b(:,1))
hold on
scatter3(my,mx,realD(:,1))


[XmiS0a XmiS1a XmiS2a XmiS3a, W0, W1, W2, W3]=multiloc_030117(   [ cell2mat(  c5unil(i) ), cell2mat(  c2unil(i) ) ],cell2mat(  c3unil(i) ) , x,  [0.3 ,0.8],3,'Epan')


[W,fit]=multiloc(   [ cell2mat(  c5unil(i) ), cell2mat(  c2unil(i) ) ],cell2mat(  c3unil(i) ) , [cgridx,cgridy], [0.3 ,0.8],3,'Epan');

[W,fit]=multiloc(   [ cell2mat(  c5unil(i) ), cell2mat(  c2unil(i) ) ],cell2mat(  c3unil(i) ) , [mx,my], [0.3 ,0.8],3,'Epan');

scatter3(x(:,1),x(:,2),XmiS2a)
hold on
scatter3(mx,my,2*fit(:,5))

while(j<2)
    

%Tmon=10;
%Tmat=10;
%[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simu2(N, 1/12, 0.02, 0.03, 0.18, 0.02, j+1, 20 ) ;    
%j=160
[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simu2(N, 1/12, 0.02, 0.03, 0.18, 0.1, j*N*T+1,T, 256 ) ;            

%OUTPUTS
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

%INPTUS
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
    
%i=1
%[XmiS0m XmiS1m XmiS2m XmiS3m]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,0.1,0.2,1  );

%%PERFORM FPCA1 ( \tilde{M}^{(0)} )only applicable if common random grid is present
%[hX2, V2, loads,Meansmo2b,D]        =fpca1(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my)
%%PERFORM FPCA1a ( \hat{M}^{(0)} )

%%Estimate the variances once to save time
sigma=zeros(N,1);
parfor(i=1:N)
sigma(i)=varaince(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
end

[hX2b, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma,0);
[hX2a, V2a, loadsa,Meansmo2ba,Da]= fpca1_gpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,2,'cpu',sigma,0);
[hX2i]=individual(Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma);


%[hX2b, V2b, loadsb,Meansmo2bb,Db,mxb,myb]=fpca1_gpu(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,0,0,1,128);
%[hX2a, V2a, loadsa,Meansmo2ba,Da,mxa,mya]= fpca1_gpu(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,0,0,2,128);


%[hX2a, V2a, loadsa,Meansmo2ba,Da]        =fpca1b(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my)
%%PERFORM FPCA2  ( \hat{M}^{(d)} )


%[hX2a, V2a, loadsa,Meansmo2ba,Da]   =fpca2(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my)

%%OUTPUTS
%Low dimensional decomposition using L Dimensions
%Eigenvectors of second derivative
%Corresponding loadings
%Mean curve
%Eigenvalues

%hXm1=[hXm1, mean( mean( (hX2-realD).^2) )]
hXm1b=[hXm1b,  mean( (hX2b-realD).^2) ];
hXm1a=[hXm1a,  mean( (hX2a-realD).^2) ];
hXm1i=[hXm1i,  mean( (hX2i-realD).^2) ];
j=j+1
end
 
scatter3(my,mx,hX2b(:,2))
hold on
scatter3(my,mx,hX2a(:,2))
hold on
scatter3(my,mx,hX2i(:,2))
hold on
scatter3(my,mx,realD(:,2))
% 
% hold on
% for (i=1:N)
% scatter3(my,mx,realD(:,i))
% end

mean(hXm1b)   %% M_0 Method
mean(hXm1a)   %% M_d Method
mean(hXm1i)   %% Individual

var(hXm1b)   %% M_0 Method
var(hXm1a)   %% M_d Method
var(hXm1i)   %% Individual

median(hXm1b)   %% M_0 Method
median(hXm1a)   %% M_d Method
median(hXm1i)   %% Individual

iqr(hXm1b)   %% M_0 Method
iqr(hXm1a)   %% M_d Method
iqr(hXm1i)   %% Individual


meanM=[meanM , [mean(hXm1b) mean(hXm1a) mean(hXm1i)]'  ]
varV= [varV , [var(hXm1b) var(hXm1a) var(hXm1i)]' ]  %% Individual
medM= [medM , [median(hXm1b) median(hXm1a) median(hXm1i)]' ]
iqrM= [iqrM , [iqr(hXm1b) iqr(hXm1a) iqr(hXm1i)]' ]

end

end

%%true Error
%scatter3(cell2mat(c5unil(1)),cell2mat(c2unil(1)),cell2mat(c3unil(1)))
%hold on
%scatter3(cell2mat(c5unil(1)),cell2mat(c2unil(1)),cell2mat(c3unilr(1)))


% % % N=50
% % % Tmon=50
% % % Tmat=6
% % % [Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simu(N, 1/12, 0.02, 0.03, 0.18, 0.02, j+1, Tmon,Tmat)             
% % % 
% % % tic
% % % [hX2b, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my,1,0,'gpu');
% % % toc
% % % 
% % % tic
% % % [hX2b, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my,1,0,'gpa');
% % % toc
% % % 
% % % tic
% % % [hX2bb, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my,1,0,'cpu');
% % % toc





[hX2b, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma,0);
