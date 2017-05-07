
%  Example file for how to use FPCder.m
%  Random design, 200 subjects with 1 to 8 measurements on each subject.
%  Time interval is [0,1]. The number of measurements for each subject is
%  uniformly distributed on {1,2,...,8}. The timepoints for each subject
%  are distributed as Beta(0.4,0.3).
%
%  In this example, the goal is to predict trajectories and their first two derivatives.

%generate data set
clear all;

addpath(genpath('PACE_matlab-master\release2.17\PACE')); 
addpath(genpath('PACE_matlab-master\release2.17\PACE\PACE-DER')); 


p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
rand('twister',sum(10000*clock));
mtp = 8;        %at most 8 repeated measurements in the simulated data
ncohort=50;    %200 subjects in the simulated data
lint=1;
y=cell(1,ncohort);
t=cell(1,ncohort);
newy = y;
newt = t;
xi=zeros(ncohort,2);
ngrid = 100;

%Case iii) regular data with missing values (regular = 1)   
%ni_missing = zeros(1,ncohort);
%for i = 1:ncohort
%  ni_missing(i) = poissrnd(mysample(1:mtp,1,0),1,1); 
%  if ni_missing(i) >= mtp
%     ni_missing(i) = mtp-1;
% end
%end


%%%Define Basisfuntions (Legrende)
syms x % define symbolic math variables
% g1p(x)=0.5*(3*x.^2-1)
% g2p(x)=1/8*(63*x.^5-70*x.^3+15*x)
% 
% g1(x)=int(g1p(x)/sqrt( int(g1p(x)^2, x,-1,1))  )
% g2(x)=int(g2p(x)/sqrt( int(g2p(x)^2, x,-1,1))  )

%g1(x)=int(0.5*(3*x.^2-1)  )
%g2(x)=int(1/8*(63*x.^5-70*x.^3+15*x) )

 g1(x)=sqrt(2)*int(0.5*(3*x.^2-1)/sqrt( int(0.5^2*(3*x.^2-1)^2, x,-1,1))  )
 g2(x)=sqrt(2)*int(1/8*(63*x.^5-70*x.^3+15*x)/sqrt( int(1/8^2*(63*x.^5-70*x.^3+15*x)^2, x,-1,1))  )

ntp=150;

 
 reps=1
 results=zeros(reps, 7)
 

 for k=1:reps
 
 
parfor i=1:ncohort

   %Case i) Sparse and irregular case (regular = 0) 
   
   %ntp=ceil(mtp*rand(1));
  
   
   t{i}=lint*rand(1,ntp);                 
   %newt{i} = lint*rand(1,ntp);
 
   %Case ii) complete balance case (regular = 2) 
   %t{i} = linspace(0,lint,mtp);                      
   %newt{i} = lint*rand(1,ntp);

   %Case iii) regular data with missing values (regular = 1)   
   %t{i} = linspace(0,lint,mtp);                    
   %newt{i} = linspace(0,lint,mtp);
   %if ni_missing(i) > 0
   %  id_missing = mysample(1:mtp, ni_missing(i),0);
   %  t{i}(id_missing) = [];
   %  newt{i}(mysample(1:mtp,ni_missing(i),0)) = [];
   %end

   xi(i,:)=[3*randn(1) 2*randn(1)];     %generate 2 Principal components
                                        %1st PC score: N(0,9), 2nd PC score: N(0,4)
   %generate the repeated measurements with measurement errors
   
   %y{i}=genMeanFun_1(t{i},0)+xi(i,:)*genEigenFun_1(t{i},0)+randn(1,length(t{i}));
   %newy{i} = genMeanFun_1(newt{i},0) + xi(i,:)*genEigenFun_1(newt{i},0)+randn(1,length(newt{i}));
   
   basis=[g2(2*(t{i}-0.5)); g1(2*(t{i}-0.5))];
   %plot(t{i},basis)
   y{i}=double( xi(i,:)*basis+0.1*randn(1,length(t{i})));

   %measurement error is distributed as N(0,1)
end


%%Perform Müller
%p = setDerOptions('yname','x','selection_k', 'FVE','FVE_threshold', 0.85,'screePlot',1, 'corrPlot',1,...
%		  'designPlot',1,'numBins',0, 'nder',0:2,'newdata',linspace(0,1,ngrid), 'ngrid',ngrid,'verbose','on');  

      
p = setDerOptions('selection_k', 2,'screePlot',0, 'corrPlot',0,...
         'designPlot',0,'numBins',128, 'kernel', 'epan',  'nder',0:2, 'verbose','off')     

% p = setDerOptions('yname','x','selection_k', 2,'FVE_threshold', 0.85,'screePlot',0, 'corrPlot',0,...
% 		  'designPlot',0,'numBins',128, 'kernel', 'epan', 'nder',0:1,'newdata',linspace(0,1,ngrid), 'ngrid',ngrid,'verbose','off');  

[yy] = FPCder(y,t,p);
ypred = getVal(yy,'y_pred');
y0 = ypred{1};     %estimated curve
y1 = ypred{2};     %estimated first derivative
%y2 = ypred{3};     %estimated second derivative
out1 = getVal(yy,'out1');
%out1 = (1:100)/100


%xtrue = repmat(genMeanFun_1(out1,0),ncohort,1)+xi*genEigenFun_1(out1,0);    %matrix of true x(t)
%xtrue1 = repmat(genMeanFun_1(out1,1),ncohort,1)+xi*genEigenFun_1(out1,1); %matrix of true x'(t)
%xtrue2 = repmat(genMeanFun_1(out1,2),ncohort,1)+xi*genEigenFun_1(out1,2); %matrix of true x''(t)

basis=[subs(g2,x,2*(out1-0.5)); subs(g1,x,2*(out1-0.5))];
basisd1=[subs(diff(g2,1),x,2*(out1-0.5)); subs(diff(g1,1),x,2*(out1-0.5))];
basisd2=[subs(diff(g2,2),x,2*(out1-0.5)); subs(diff(g1,2),x,2*(out1-0.5))];

xtrue  = xi*double( basis   );    %matrix of true x(t)
xtrue1 = xi*double( basisd1 ); %matrix of true x'(t)
xtrue2 = xi*double( basisd2 ); %matrix of true x''(t)



addpath(genpath('base')); 
Yc=[]
Xc=[]
N=length(y)
for i=1:N
    Xc{i}= cell2mat(t(i))'; 
    Yc{i}= cell2mat(y(i))';
end


%%Univariate estimate for comparison      
%%Derive Bandwidth
order=3;

    sigma=zeros(N,1);
    parfor i=1:N
        sigma(i)=varaince(cell2mat( Xc(i) ),cell2mat(  Yc(i) ) );
    end
    
    sigma=repmat(1,N,1)*mean(sigma)

    Lpx=zeros(size(out1',1),1);
    for j=1:(order+2)
        Lpx      =[Lpx, factorial(j)*out1'.^max(0,j-order )*min(1, max(0,j-order+1 ) )  ]; 
    end
    
    scores=[];
    for i=1:N
        xachs   =cell2mat(  Xc(i) );
        Lp=ones(size(xachs,1),1);
        for j=1:(order+2)
            Lp      =[Lp, xachs.^j];
        end
        Vxre    =cell2mat(  Yc(i) );
        scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
    end
    regsx=Lpx*scores;
g=1
    ha=zeros(g,1);
        ha= mean(sigma)./( ntp.*( mean( (regsx).^2  )));

        
%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE
order=2
der=1
pos2=der+1
int=((-300:300)/100)';
eqkernel=equivkernel(int,order,'Gauss');
Cdp=factorial(order+1)^2*(2*sum(der)+1)*range(int)*mean(eqkernel(pos2,:).^2)/(2*(order+1-sum(der) )*range(int)^2*mean(int'.^(order+1).*eqkernel(pos2,:) )^2 ) ;
Cdp^(1/(2*order+3))
Hd=diag( mean(Cdp*ha').^(1/(2*order+3) ) )

fit_ind2=[]
for i=1:N
[W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Yc(i) ),out1' ,Hd,order,'Gauss'  );
fit_ind2=[fit_ind2, fit(:,2) ];
end

order=1
der=0
pos2=der+1
int=((-300:300)/100)';
eqkernel=equivkernel(int,order,'Gauss');
Cdp=factorial(order+1)^2*(2*sum(der)+1)*range(int)*mean(eqkernel(pos2,:).^2)/(2*(order+1-sum(der) )*range(int)^2*mean(int'.^(order+1).*eqkernel(pos2,:) )^2 ) ;
Cdp^(1/(2*order+3))
H0=diag( mean(Cdp*ha').^(1/(2*order+3) ) )
fit_ind=[]
for i=1:N
[W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Yc(i) ),out1' ,H0,order,'Gauss'  );
fit_ind=[fit_ind, fit(:,1) ];
end
     
[hX2b, V2b, loadsb,Meansmo2bb,Db2,xo]=fpca(2,1,'M_0',Xc,Yc,out1',sigma,0,0,1);
[hX2bd, V2bd, loadsb,Meansmo2bb,Db2,xo]=fpca(2,1,'M_d',Xc,Yc,out1',sigma,0,0,1);
[hX1b, V1b, loadsb,Meansmo2bb,Db1,xo]=fpca(2,0,'M_0',Xc,Yc,out1',sigma,0,0,1);
 


results(k,:)   = [mean( mean( (xtrue1-hX2b').^2 ) ) ,mean( mean( (xtrue1-hX2bd').^2 ) ), mean( mean( (xtrue1-fit_ind2').^2 ) ), mean( mean( (xtrue1-cell2mat(y1')).^2 ) ) , mean( mean( (xtrue-hX1b').^2 ) ) , mean( mean( (xtrue-fit_ind').^2 ) ) ,mean( mean( (xtrue-cell2mat(y0')).^2 ) ) ];

%%%Define Basisfuntions (Legrende)
%syms x y% define symbolic math variables
%g1(x)= x+ sin y
%g2(x)= x+ cos y


%[hX2b, V2b, loadsb,Meansmo2bb,Db,x_out]=fpca(2,2,'M_0',Xc,Yc,out1',1,0,0,1);
 
 end
 
dataset({results 'New M0 Method d=1' ,'New Md Method d=1',  'Individual d=1', 'Mueller Method d=1', 'New Method d=0', 'Individual d=0', 'Mueller Method d=0'})


i=5

hold off
plot(out1,xtrue1(i,:),out1,hX2b(:,i),out1,fit_ind2(:,i),out1,y1{i},'b-');
hold on
title(['Subject ' num2str(i)]);
xlabel('t');
ylabel('X(t)');
legend('true','pred','pred indiv','pred Müller','Location','Best');



phi = getVal(yy,'phi');        %estimated eigenfunctions
phi0 = phi{1};                 %estimated eigenfunctions
phi1 = phi{2};                 %estimated 1st derivative of eigenfunctions


hold off
plot(out1,xtrue1(i,:),out1,hX2b(:,i),out1,fit_ind2(:,i),out1,y1{i},'b-');
hold on
title(['Subject ' num2str(i)]);
xlabel('t');
ylabel('X(t)');
legend('true','pred','pred indiv','pred Müller','Location','Best');

basisc=double(basis)'


plot(out1,V1b(:,2))
hold on
plot(out1,phi{1}(:,2))
hold on
plot(out1,basisc(:,2))


plot(out1,V1b(:,2))
hold on
plot(out1,phi{1}(:,2))
hold on
plot(out1,basisc(:,2)/sqrt( mean(basisc(:,1).^2) ))

plot(out1,-V1b(:,1))
hold on
plot(out1,phi{1}(:,1))
hold on
plot(out1,basisc(:,1)/sqrt( mean(basisc(:,2).^2) ) )



%%%Second eigenfunction

basiscd1=double(basisd1)'
plot(out1, -V2bd(:,2) )
hold on
plot(out1,basiscd1(:,2))

plot(out1,-V2bd(:,1) )
hold on
plot(out1,basiscd1(:,1))



plot(out1,phi{2}(:,2))
hold on
plot(out1,phi{2}(:,1) )
hold on
phi0 = phi{1};                 %estimated eigenfunctions
phi1 = phi{2};                 %estimated 1st derivative of eigenfunctions
phi2 = phi{3};                 %estimated 2nd derivative of eigenfunctions

