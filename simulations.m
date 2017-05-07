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
medianPhi=[ ]


m=0
for (N=[10 25])
for (T=[50 250])
%for (N=[10])
%for (T=[50])

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

j=1;

while(j<501)
    

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
Xc=[]
Yc=[]
for i=1:N
    Xc{i}= [ cell2mat(c2unil(i)),cell2mat(c5unil(i)) ]
    Yc{i}= cell2mat(c3unil(i));
end


    sigma=zeros(N,1);
    parfor i=1:N
        sigma(i)=varaince(cell2mat( Xc(i) ),cell2mat(  Yc(i) ) );
    end

   
[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(2,[2 0],'M_0',Xc,Yc,[mx,my],sigma,0,[min(x1minc),min(x2minc)],[max(x1maxc),max(x2maxc)],'Gauss','Gauss');
[hX2a, V2a, loadsa,Meansmo2ba,Da,x]=fpca(2,[2 0],'M_d',Xc,Yc,[mx,my],sigma,0,[min(x1minc),min(x2minc)],[max(x1maxc),max(x2maxc)],'Gauss','Gauss');
[hX2i]=individual(Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c2unil,c5unil,c3unil,c3unilr,N,mx,my,1,'cpu',sigma,3,'Gauss');


hXm1b=[hXm1b,  mean( (hX2b-realD).^2)./(mean(mean(realD.^2))) ];
hXm1a=[hXm1a,  mean( (hX2a-realD).^2)./(mean(mean(realD.^2))) ];
hXm1i=[hXm1i,  mean( (hX2i-realD).^2)./(mean(mean(realD.^2))) ];


%%Compare Eigenfunctions
%Get eigenvalues
X0mean=cell2mat(c3unilr) -repmat( mean(cell2mat(c3unilr),2) ,1,N);
Xdmean=realD -repmat( mean(realD,2) ,1,N);

[VrealM0, DrealM0]=eig(  prod(range([mx,my]))/length(mx)*X0mean'*X0mean  );
[DrealM0,I]=sort(diag(DrealM0),'descend');
DrealM0=diag(DrealM0);
VrealM0=VrealM0(:,I);
pc0Real=Xdmean * VrealM0(:,1:2) * diag( diag(DrealM0(1:2,1:2)).^(-1/2) );

[VrealMd, DrealMd]=eig(  prod(range([mx,my]))/length(mx)*Xdmean'*Xdmean  );
[DrealMd,Id]=sort(diag(DrealMd),'descend');
DrealMd=diag(DrealMd);

VrealMd=VrealMd(:,Id);
pcdReal=Xdmean * VrealMd(:,1:2) * diag( diag(DrealMd(1:2,1:2)).^(-1/2) );

% 
% Ph1ia=[Ph1ia , mse_ord(pcdReal(:,1),V2a(:,1),33)./mean(pcdReal(:,1).^2) ]
% Ph2ia=[Ph2ia , mse_ord(pcdReal(:,2),V2a(:,2),33)./mean(pcdReal(:,2).^2) ]
% Ph1ib=[Ph1ib , mse_ord(pc0Real(:,1),V2b(:,1),33)./mean(pc0Real(:,1).^2) ]
% Ph2ib=[Ph2ib , mse_ord(pc0Real(:,2),V2b(:,2),33)./mean(pc0Real(:,2).^2) ]


Ph1ia=[Ph1ia , mse_ord(pcdReal(:,1)/sqrt(mean(pcdReal(:,1).^2))  ,V2a(:,1)/sqrt(mean(V2a(:,1).^2)),33)./mean(pcdReal(:,1).^2)  ]
Ph2ia=[Ph2ia , mse_ord(pcdReal(:,2)/sqrt(mean(pcdReal(:,2).^2))  ,V2a(:,2)/sqrt(mean(V2a(:,2).^2)),33)./mean(pcdReal(:,2).^2) ]
Ph1ib=[Ph1ib , mse_ord(pc0Real(:,1)/sqrt(mean(pc0Real(:,1).^2)),V2b(:,1)/  sqrt(mean(V2b(:,1).^2)),33)./mean(pc0Real(:,1).^2) ]
Ph2ib=[Ph2ib , mse_ord(pc0Real(:,2)/sqrt(mean(pc0Real(:,2).^2)),V2b(:,2)/  sqrt(mean(V2b(:,2).^2)),33)./mean(pc0Real(:,2).^2) ]

prod(range(x))*mean( pcdReal(:,1).^2)
prod(range(x))*mean( pc0Real(:,1).^2 )
prod(range(x))*mean( V2b(:,1).^2 )

j=j+1
end
 
FFun1= TriScatteredInterp(mx,my,hX2b(:,1), 'linear')
% define axis
axis_x=linspace(min(mx),max(mx),50); % maturity
axis_y=linspace(min(my),max(my),50); % moneyness
[Xa,Ya] = meshgrid(axis_x,axis_y);
Za1b=FFun1(Xa,Ya);
    surfl(Xa,Ya,Za1b);
    set(gca,'XTick',[0.50 0.75 1 1.25 1.50])
    set(gca,'YTick',[0 0.25 0.50 0.75 1 ])
    xh=xlabel('Moneyness')
    yh=ylabel('Maturity')
    zh=zlabel('$$\hat{ \gamma}^{(d)}_{2,T}$$')
    set(zh,'Interpreter','latex')
    zlim([-120,100])
    set([xh,yh,zh],'fontsize',10)
    set(gca, 'fontsize',10)
    
    
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


meanM=[meanM , [T N mean(hXm1b) mean(hXm1a) mean(hXm1i)]'  ]
varV= [varV , [T N var(hXm1b) var(hXm1a) var(hXm1i)]' ]  %% Individual
medM= [medM , [T N median(hXm1b) median(hXm1a) median(hXm1i)]' ]
iqrM= [iqrM , [T N iqr(hXm1b) iqr(hXm1a) iqr(hXm1i)]' ]

meanPhi=[meanPhi ,[T N mean(Ph1ib) mean(Ph2ib) mean(Ph1ia) mean(Ph2ia) ]' ]
medianPhi=[medianPhi ,[T N median(Ph1ib) median(Ph2ib) median(Ph1ia) median(Ph2ia) ]' ]


end

end


scatter3(my,mx,hX2b(:,1),'filled')
hold on
scatter3(my,mx,hX2a(:,1),'filled')
hold on
scatter3(my,mx,hX2i(:,1),'filled')
hold on
scatter3(my,mx,realD(:,1))



%scatter3(my,mx,1/sqrt(Db(2))*V2b(:,2),'filled')

pcdReal(:,1)/mean(pcdReal(:,1).^2)  ,V2a(:,1)/mean(V2a(:,1).^2)

pc0Real(:,1)/mean(pc0Real(:,1).^2),V2b(:,1)/mean(V2b(:,1).^2)


 scatter3(my,mx,V2b(:,1)/sqrt( mean(V2b(:,1).^2)) ,'filled')
 hold on
scatter3(my,mx, pc0Real(:,1)/sqrt( mean(pc0Real(:,1).^2)))


Db(1:5)
diag( DrealMd(1:5,1:5) )

csvwrite('meang2.csv', meanM')
csvwrite('varVg2.csv', varV')
csvwrite('medMg2.csv', medM')
csvwrite('iqrMg2.csv', iqrM')


dlmwrite('meanPhi2.csv', [ meanPhi(1:2,:) ; round( meanPhi(3,:)*mean( pc0Real(:,1).^2),2); round( meanPhi(4,:)*mean( pc0Real(:,2).^2),2) ; round( meanPhi(5:6,:),2)]', 'precision','%.4f')
dlmwrite('meedianPhi2.csv', [ medianPhi(1:2,:) ; round( medianPhi(3,:)*mean( pc0Real(:,1).^2),2 ); round( medianPhi(4,:)*mean( pc0Real(:,2).^2),2 ) ; round( medianPhi(5:6,:),2) ]', 'precision','%.4f')

scatter3(my,mx, Meansmo2bb )
scatter3(my,mx, Meansmo2ba )