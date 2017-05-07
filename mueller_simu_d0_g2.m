% % % clear all;
% % % Derive optimal bw for M
% % % syms t
% % % syms C1 %d!/(p+1)!
% % % syms d
% % % syms p
% % % syms C2 %int int (u^(p+1) Kd^2)X^(p+1)X dz
% % % syms C3 %2 d!^4 int( sigma^4 2*sigma^2 X(x)^2) f(x)^2 dx int K^C^2 du
% % % syms T
% % % syms b
% % % 
% % % MSE= (2*C1*b^(p+1-2)*C2)^2 + b^(-4*d-1)*T^(-2)*C3
% % % 
% % % diff(MSE,b)
% % % 
% % % b_stern=solve(diff(MSE,b),b)
% % % 
% % % 
% % % g=C1*b+C2+C3*3
% % % x=solve(g,b) 

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

addpath(genpath('base')); 

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
t2=cell(1,ncohort);
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
syms x y % define symbolic math variables

g1(x,y)=sin( x*(4*pi) )*sin( y*(4*pi) )
g2(x,y)=cos(2*x*pi)*cos(2*y*pi)


ntp=150;


 reps=1
 results=zeros(reps, 7)
 
 resultsEV=zeros(reps, 8)
 

 for k=1:reps
 
 
parfor i=1:ncohort

   
   t{i}=[lint*rand(1,ntp)' lint*rand(1,ntp)'];

   
   xi(i,:)=[3*randn(1) 2*randn(1)];     %generate 2 Principal components
      
   %basis=[g2(2*(t{i}-0.5)); g1(2*(t{i}-0.5))];
   basis=[g2((t{i}')); g1((t{i}))];
   %plot(t{i},basis)
   y{i}= xi(i,:)*double(basis)+2*randn(1,length(t{i}));

   %measurement error is distributed as N(0,1)
end


%%Perform Müller
%p = setDerOptions('yname','x','selection_k', 'FVE','FVE_threshold', 0.85,'screePlot',1, 'corrPlot',1,...
%		  'designPlot',1,'numBins',0, 'nder',0:2,'newdata',linspace(0,1,ngrid), 'ngrid',ngrid,'verbose','on');  

      
p = setDerOptions('selection_k', 2,'screePlot',0, 'corrPlot',0,...
         'designPlot',0,'numBins',128, 'kernel', 'epan',  'nder',0:2, 'verbose','off')     

% p = setDerOptions('yname','x','selection_k', 2,'FVE_threshold', 0.85,'screePlot',0, 'corrPlot',0,...
% 		  'designPlot',0,'numBins',128, 'kernel', 'epan', 'nder',0:1,'newdata',linspace(0,1,ngrid), 'ngrid',ngrid,'verbose','off');  

%[yy] = FPCder(y,t,p);
%ypred = getVal(yy,'y_pred');
%y0 = ypred{1};     %estimated curve
%y1 = ypred{2};     %estimated first derivative
%y2 = ypred{3};     %estimated second derivative
%out1 = getVal(yy,'out1');
%out1 = (1:100)/100


%xtrue = repmat(genMeanFun_1(out1,0),ncohort,1)+xi*genEigenFun_1(out1,0);    %matrix of true x(t)
%xtrue1 = repmat(genMeanFun_1(out1,1),ncohort,1)+xi*genEigenFun_1(out1,1); %matrix of true x'(t)
%xtrue2 = repmat(genMeanFun_1(out1,2),ncohort,1)+xi*genEigenFun_1(out1,2); %matrix of true x''(t)

% basis=[subs(g2,x,2*(out1-0.5)); subs(g1,x,2*(out1-0.5))];
% basisd1=[subs(diff(g2,1),x,2*(out1-0.5)); subs(diff(g1,1),x,2*(out1-0.5))];
% basisd2=[subs(diff(g2,2),x,2*(out1-0.5)); subs(diff(g1,2),x,2*(out1-0.5))];

basis=[subs(g2,x,(out1)); subs(g1,x,(out1))];
basisd1=[subs(diff(g2,1),x,(out1)); subs(diff(g1,1),x,(out1))];
basisd2=[subs(diff(g2,2),x,(out1)); subs(diff(g1,2),x,(out1))];

xtrue  = xi*double( basis   );    %matrix of true x(t)
xtrue1 = xi*double( basisd1 ); %matrix of true x'(t)
xtrue2 = xi*double( basisd2 ); %matrix of true x''(t)



addpath(genpath('base')); 
Yc=[]
Xc=[]
N=length(y)
for i=1:N
    Xc{i}= [ cell2mat(t(i))', cell2mat(t2(i))']; 
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
[hX2bd, V2bd, loadsb,Meansmo2bb,Db2d,xo]=fpca(2,1,'M_d',Xc,Yc,out1',sigma,0,0,1);
[hX1b, V1b, loadsb,Meansmo2bb,Db1,xo]=fpca(2,0,'M_0',Xc,Yc,out1',sigma,0,0,1);
 


results(k,:)   = [mean( mean( (xtrue1-hX2b').^2 ) ) ,mean( mean( (xtrue1-hX2bd').^2 ) ), mean( mean( (xtrue1-fit_ind2').^2 ) ), mean( mean( (xtrue1-cell2mat(y1')).^2 ) ) , mean( mean( (xtrue-hX1b').^2 ) ) , mean( mean( (xtrue-fit_ind').^2 ) ) ,mean( mean( (xtrue-cell2mat(y0')).^2 ) ) ];


basisc=double(basis)'
basiscd1=double(basisd1)'
phi = getVal(yy,'phi');        %estimated eigenfunctions

resultsEV(k,:)   = [mse_ord(basisc(:,1),V1b(:,1),33), mse_ord(basisc(:,2),V1b(:,2),33), mse_ord(basisc(:,1),phi{1}(:,1),33) ,mse_ord(basisc(:,2),phi{1}(:,2),33),mse_ord(basiscd1(:,1),V2b(:,1),33),mse_ord(basiscd1(:,2),V2b(:,2),33),mse_ord(basiscd1(:,1),phi{2}(:,1),33),mse_ord(basiscd1(:,2),phi{2}(:,2),33)]

%%%Define Basisfuntions (Legrende)
%syms x y% define symbolic math variables
%g1(x)= x+ sin y
%g2(x)= x+ cos y


%[hX2b, V2b, loadsb,Meansmo2bb,Db,x_out]=fpca(2,2,'M_0',Xc,Yc,out1',1,0,0,1);
 
 end
 
hX2b 
V2b 
Meansmo2bb

dataset({results 'New M0 Method d=1' ,'New Md Method d=1',  'Individual d=1', 'Mueller Method d=1', 'New Method d=0', 'Individual d=0', 'Mueller Method d=0'})

dataset({resultsEV 'New M0 Method phi1' ,'New M0 Method phi2', 'Mueller phi1', 'Mueller phi2', 'New M0 Method d_phi1', 'New M0 Method d_phi2', 'Mueller d_phi1' , 'Mueller d_phi2'})


cell2mat( Xc )


i=2

hold off
plot(out1,xtrue(i,:),out1,hX1b(:,i),out1,fit_ind(:,i),out1,y0{i},'b-') ;
hold on
scatter(cell2mat(Xc(i)),cell2mat(Yc(i)))
title(['Subject ' num2str(i)]);
xlabel('t');
ylabel('X(t)');
legend('true','pred','pred indiv','pred Müller','observed','Best');



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
plot(out1,-V1b(:,1))
hold on
plot(out1,phi{1}(:,1))
hold on
plot(out1,-basisc(:,1))


plot(out1,V1b(:,2))
hold on
plot(out1,phi{1}(:,2))
hold on
plot(out1,basisc(:,2))

%%%Second derivative of eigenfunction 

basiscd1=double(basisd1)'
plot(out1, -V2b(:,1) )
hold on
plot(out1,-phi{2}(:,1))
hold on
plot(out1,basiscd1(:,1))

plot(out1,V2b(:,2) )
hold on
plot(out1,phi{2}(:,2))
hold on
plot(out1,basiscd1(:,2))


% 
% %%%Second eigenfunction 
% 
% basiscd1=double(basisd1)'
% plot(out1, -V2bd(:,1) )
% hold on
% plot(out1,basiscd1(:,2)/(2*pi ) )
% 
% plot(out1,-V2bd(:,2) )
% hold on
% plot(out1,basiscd1(:,1)/(4*pi) )
% 


%TESTAREA IS POOLING THE SAME ???

H0=diag( mean(1/N*Cdp*ha').^(1/(2*order+3) ) )


fit_ind=[]
for i=1:N
[W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Yc(i) ),out1' ,H0,order,'Gauss'  );
fit_ind=[fit_ind, fit(:,1) ];
end


plot(out1, mean(fit_ind,2) )
hold on
%%compare with pooled mean

H0N=diag( mean(Cdp*ha').^(1/(2*order+3) ) )

[W_dm,fitpool]=multiloc( cell2mat( Xc' ) , cell2mat(  Yc' ),out1' ,H0/1.1,order,'Gauss'  );
fitpool(:,1)

plot(out1,fitpool(:,1)')




scatter(cell2mat( Xc' ), cell2mat(  Yc' ))
hold on 
plot(out1, fit_ind )
