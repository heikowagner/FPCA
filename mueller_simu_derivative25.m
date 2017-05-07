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

addpath(genpath('simulation')); 
addpath(genpath('base')); 
%addpath(genpath('I:\OwnCloud2\Uni\Programs\Quantanet Version_neu')); 
%addpath(genpath('I:\OwnCloud2\Uni\Programs\Quantanet Version_GPU\')); 
%addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_GPU'));  
%addpath(genpath('C:\Users\Heiko Wagner\Desktop\Quantanet Version_GPU')); 
addpath(genpath('PACE_matlab-master\release2.17\PACE')); 
addpath(genpath('PACE_matlab-master\release2.17\PACE\PACE-DER')); 


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
meanPhi=[]
varPhi=[]

Scores=[]
eigenvalues=[]
Medeigenvalues=[]

medianPhi=[ ]

meanScores=[  ]
medianScoresa=[  ]


m=0
%for (N=[10 25 200])
%for (T=[15 100 250])
for (N=[25 50 200])
for (T=[25 50 200])

m=m+1    

hXm1b=[];
hXm1a=[];
hXm1i=[];
Ph1ia=[];
Ph2ia=[];
Ph1ib=[];
Ph2ib=[];
Ph1im=[];
Ph2im=[];

eigenA= [ ];
eigenB= [ ];
eigenMU=[ ];

Scoresa=[ ];
Scoresb=[ ];
Scoresmu=[ ];

j=1;

while(j<501)
%while(j<3)
%while(j<10)
    

%Tmon=10;
%Tmat=10;
%[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simu2(N, 1/12, 0.02, 0.03, 0.18, 0.02, j+1, 20 ) ;    
%j=160

%[x1minc,x1maxc,cgridx,c3unil,c3unilr,realD,mx,Y_Dtrue]=simu1(N, 1/12, 0.02, 0.03, 0.05, j*N*T+1,T, 256 ) ;            
[x1minc,x1maxc,cgridx,c3unil,c3unilr,Ytrue, realD,mx,Y_Dtrue]=simu1(N, 1/12, 0.02, 0.03, 0.005, 100*j*N*T+1,T, 256 ) ;            

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
Xc2=[]
Yc2=[]
for i=1:N
    Xc2{i}=  [ cgridx(:,i)' ];
    Yc2{i}=  [ c3unil(:,i)' ];
end


    
% p = setDerOptions('selection_k', 2,'screePlot',0, 'corrPlot',0,...
%             'designPlot',0,'numBins',15, 'kernel', 'epan',  'nder',0:2, 'verbose','off')     

p = setDerOptions('selection_k', 2,'screePlot',0, 'corrPlot',0,...
            'designPlot',0,'numBins',15, 'nder',0:2, 'verbose','off')     

      
[yy] = FPCder(Yc2,Xc2,p);
ypred = getVal(yy,'y_pred');
y0 = ypred{1};     %estimated curve
y1 = ypred{2};     %estimated first derivative
y2 = ypred{3};     %estimated second derivative
out1 = getVal(yy,'out1');
phi = getVal(yy,'phi');        %estimated eigenfunctions
phi0 = phi{1};                 %estimated eigenfunctions
phi1 = phi{2};                 %estimated 1st derivative of eigenfunctions
phi2 = phi{3};                 %estimated 2st derivative of eigenfunctions
evMU=getVal(yy,'lambda')
    
Xc=[]
Yc=[]
for i=1:N
    Xc{i}=  [ cgridx(:,i) ];
    Yc{i}=  [ c3unil(:,i) ];
end

    sigma=zeros(N,1);
    parfor i=1:N
        sigma(i)=varaince(cell2mat( Xc(i) ),cell2mat(  Yc(i) ) );
    end

[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(2,2,'M_0',Xc,Yc,mx',sigma,0,min(x1minc),max(x1maxc),'Gauss','Gauss');
[hX2a, V2a, loadsa,Meansmo2ba,Da,x]=fpca(2,2,'M_d',Xc,Yc,mx',sigma,0,min(x1minc),max(x1maxc),'Gauss','Gauss');


hX2lm=[]
for i=1:N
        hX2lm=[hX2lm interp1(out1',cell2mat( y2(i) ),mx','linear','extrap')];
end

phi22=[interp1(out1',phi2(:,1),mx','linear','extrap') interp1(out1',phi2(:,2),mx','linear','extrap')];

hXm1b=[hXm1b,  mean( (hX2b-realD).^2)./mean( mean(realD.^2)) ];
hXm1a=[hXm1a,  mean( (hX2a-realD).^2)./mean( mean(realD.^2)) ];
hXm1i=[hXm1i,  mean( (hX2lm-realD).^2)./mean( mean(realD.^2)) ];
j=j+1


%%Compare Eigenfunctions
%Get eigenvalues
X0mean=c3unilr -repmat( mean(c3unilr,2) ,1,N);
Xdmean=realD -repmat( mean(realD,2) ,1,N);

[VrealM0, DrealM0]=eig(  range(x)/length(x)*X0mean'*X0mean  );
[DrealM0,I]=sort(diag(DrealM0),'descend');
DrealM0=diag(DrealM0);
VrealM0=VrealM0(:,I);
pc0Real=Xdmean * VrealM0(:,1:2) * diag( diag(DrealM0(1:2,1:2)).^(-1/2) );

[VrealMd, DrealMd]=eig(  range(x)/length(x)*Xdmean'*Xdmean  );
[DrealMd,Id]=sort(diag(DrealMd),'descend');
DrealMd=diag(DrealMd);
VrealMd=VrealMd(:,Id);
pcdReal=Xdmean * VrealMd(:,1:2) * diag( diag(DrealMd(1:2,1:2)).^(-1/2) );

%Evaluation

eigenA= [eigenA,  (diag( DrealMd(1:2,1:2) )/N-Da(1:2)).^2 ];
eigenB= [eigenB,  (diag( DrealM0(1:2,1:2) )/N-Db(1:2)).^2 ];
eigenMU=[eigenMU,  ( (  diag( DrealM0(1:2,1:2)/N)'-evMU).^2)' ];


Ph1ia=[Ph1ia , mse_ord(pcdReal(:,1),V2a(:,1),33)./(mean(pcdReal(:,1).^2)) ]
Ph2ia=[Ph2ia , mse_ord(pcdReal(:,2),V2a(:,2),33)./(mean(pcdReal(:,2).^2)) ]
Ph1ib=[Ph1ib , mse_ord(pc0Real(:,1),V2b(:,1),33)./(mean(pc0Real(:,1).^2)) ]
Ph2ib=[Ph2ib , mse_ord(pc0Real(:,2),V2b(:,2),33)./(mean(pc0Real(:,2).^2)) ]
Ph1im=[Ph1im , mse_ord(pc0Real(:,1),phi22(:,1),33)./(mean(pc0Real(:,1).^2)) ]
Ph2im=[Ph2im , mse_ord(pc0Real(:,2),phi22(:,2),33)./(mean(pc0Real(:,2).^2)) ]
      
Scoresa=[Scoresa; mse_ord( VrealM0(:,1:2) * diag( diag(DrealM0(1:2,1:2)).^(-1/2) ) ,loadsa(1:2,:)')]
Scoresb=[Scoresb; mse_ord(VrealMd(:,1:2) * diag( diag(DrealMd(1:2,1:2)).^(-1/2) ) ,loadsa(1:2,:)')]
Scoresmu=[Scoresmu; mse_ord(VrealM0(:,1:2) * diag( diag(DrealM0(1:2,1:2)).^(-1/2) ) ,getVal(yy,'xi_est'))]
end


%plot( phi22*getVal(yy,'xi_est')'-hX2lm  +repmat( mean(hX2lm,2) ,1,200) )


eigenvalues=[eigenvalues, [T N mean(eigenB') mean( eigenMU')]' ]
Medeigenvalues=[Medeigenvalues, [T N median(eigenB') median( eigenMU')]' ]


%%Comarae Means
meanM=[meanM , [T N mean(hXm1b) mean(hXm1a) mean(hXm1i)]'  ]
varV= [varV , [T N var(hXm1b) var(hXm1a) var(hXm1i)]' ]  %% Individual
medM= [medM , [T N median(hXm1b) median(hXm1a) median(hXm1i)]' ]
iqrM= [iqrM , [T N iqr(hXm1b) iqr(hXm1a) iqr(hXm1i)]' ]

meanPhi=[meanPhi ,[T N mean(Ph1ib) mean(Ph2ib) mean(Ph1ia) mean(Ph2ia) mean(Ph1im) mean(Ph2im)]' ]
medianPhi=[medianPhi ,[T N median(Ph1ib) median(Ph2ib) median(Ph1ia) median(Ph2ia) median(Ph1im) median(Ph2im)]' ]

meanScores=[medianScoresa ,[T N mean(Scoresa) mean(Scoresb) mean(Scoresmu)]' ]
medianScoresa=[medianScoresa ,[T N median(Scoresa) median(Scoresb)  median(Scoresmu)]' ]

end

end

dataset({meanM' 'T' 'N' 'New_M0_Method_d=2' ,'New_Md_Method_d=2', 'Mueller_Method_d=2'})
dataset({varV' 'T' 'N' 'New_M0_Method_d=2' ,'New_Md_Method_d=2', 'Mueller_Method_d=2'})
dataset({medM' 'T' 'N' 'New_M0_Method_d=2' ,'New_Md_Method_d=2', 'Mueller_Method_d=2'})
dataset({iqrM' 'T' 'N' 'New_M0_Method_d=2' ,'New_Md_Method_d=2', 'Mueller_Method_d=2'})

meanPhi'
names=[{'T'},{'N'}, 'New_M0_Method_d=2', 'New_Md_Method_d=2', 'Mueller_Method_d=2']
meanOut= [names, meanM']

format bank
[ meanM(1:2,:) ; round( meanM(3:5,:)*10,2) ]'

dlmwrite('mean25.csv', [ meanM(1:2,:) ; round( meanM(3:5,:)*10,2) ]', 'precision','%.2f')
dlmwrite('varV25.csv', [ varV(1:2,:) ; round( varV(3:5,:)*10^3,2) ]', 'precision','%.2f')
dlmwrite('medM25.csv', [ medM(1:2,:) ; round( medM(3:5,:)*10^1,2) ]', 'precision','%.2f')
dlmwrite('iqrM25.csv', [ iqrM(1:2,:) ; round( iqrM(3:5,:)*10^1,2) ]', 'precision','%.2f')
dlmwrite('meanPhi25.csv', [ meanPhi(1:2,:) ; round( meanPhi(3:8,:)*10^1,2) ]', 'precision','%.2f')
dlmwrite('meedianPhi25.csv', [ medianPhi(1:2,:) ; round( medianPhi(3:8,:)*10^1,2) ]', 'precision','%.2f')
dlmwrite('eigenvalues25.csv', [ eigenvalues(1:2,:) ; round( eigenvalues(3:6,:)*10^10,2) ]', 'precision','%.2f')
dlmwrite('medianeigenvalues25.csv', [ Medeigenvalues(1:2,:) ; round( Medeigenvalues(3:6,:)*10^10,2) ]', 'precision','%.2f')
dlmwrite('meanScores25.csv', [ meanScores(1:2,:) ; round( meanScores(3:6,:)*10,2) ]', 'precision','%.2f')
dlmwrite('medianScores25.csv', [ medianScoresa(1:2,:) ; round( medianScoresa(3:6,:)*10,2) ]', 'precision','%.2f')


% csvwrite('mean.csv', [ meanM(1:2,:) ; round( meanM(3:5,:)*10,2) ]')
% csvwrite('varV.csv', [ varV(1:2,:) ; round( varV(3:5,:)*10^3,2) ]')
% csvwrite('medM.csv', [ medM(1:2,:) ; round( medM(3:5,:)*10^1,2) ]')
% csvwrite('iqrM.csv', [ iqrM(1:2,:) ; round( iqrM(3:5,:)*10^1,2) ]')
% csvwrite('meanPhi.csv', [ meanPhi(1:2,:) ; round( meanPhi(3:8,:)*10^1,2) ]')
% csvwrite('eigenvalues.csv', [ eigenvalues(1:2,:) ; round( eigenvalues(3:6,:)*10^10,2) ]')



% 
% csvwrite('mean.csv', meanM')
% csvwrite('varV.csv', varV')
% csvwrite('medM.csv', medM')
% csvwrite('iqrM.csv', iqrM')
% csvwrite('meanPhi.csv', meanPhi')
% csvwrite('eigenvalues.csv', eigenvalues')
T=25
N=25

[x1minc,x1maxc,cgridx,c3unil,c3unilr,Ytrue, realD,mx,Y_Dtrue]=simu1(N, 1/12, 0.02, 0.03, 0.005, j*N*T+1,T, 256 ) ;            

Xc=[]
Yc=[]
for i=1:N
    Xc{i}=  [ cgridx(:,i) ];
    Yc{i}=  [ c3unil(:,i) ];
end

%[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(2,2,'M_0',Xc,Yc,mx',sigma,0,min(x1minc),max(x1maxc),'Gauss','Epan');



Xc=[]
Yc=[]
for i=1:N
    Xc{i}=  [ cgridx(:,i)' ];
    Yc{i}=  [ c3unil(:,i)' ];
end

%  p = setDerOptions('selection_k', 2,'screePlot',0, 'corrPlot',0,...
%              'designPlot',0,'numBins',15,  'nder',0:2, 'verbose','off')     
% %       
%       
% [yy] = FPCder(Yc,Xc,p);
% ypred = getVal(yy,'y_pred');
% y0 = ypred{1};     %estimated curve
% y1 = ypred{2};     %estimated first derivative
% y2 = ypred{3};     %estimated second derivative
% out1 = getVal(yy,'out1');
% phi = getVal(yy,'phi');        %estimated eigenfunctions
% phi0 = phi{1};                 %estimated eigenfunctions
% phi1 = phi{2};                 %estimated 1st derivative of eigenfunctions
% phi2 = phi{3};                 %estimated 2st derivative of eigenfunctions
% evMU=getVal(yy,'lambda')
% 


%plot(out1, phi2,'Color',[0,0.1,0.9])
plot(mx, phi22,'Color',[0,0.1,0.9])
hold on
plot(mx, V2b(:,1:2),'Color','red')
hold on
plot(mx, pc0Real(:,1:2)*[1 0; 0 1 ],'Color',[0,0.8,1] )

 mse_ord(pc0Real(:,1),V2b(:,1),33) 
mse_ord(pc0Real(:,2),V2b(:,2),33) 
 mse_ord(pc0Real(:,1),phi22(:,1),70) 
mse_ord(pc0Real(:,2),phi22(:,2),70) 

