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
         %Number of Dimensions +1
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
% Author:      Heiko Wagner, 2017/01/17 
% ------------------------------------------------------------------------------ 

clear all 

%addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_neu')); 

addpath(genpath('C:\Users\Heiko\Desktop\Quantanet Version_GPU')); 
%%%READ REAL DATA
%%AMOUNT NUMBER OF DIMENSIONS USED ! (+1)
n=500;   %Grid #
m=100;
ta=1/12;   %%%Anzahl Observationen
r=0.01; %%%Riskfree Zinssatz
%mu=0.03;  %%%Drift
sigma=0.07;   %Volatility
TAU=[1/24 1/12 3/12 6/12 11/12];   %%Maturitys X2 Achse 

sigmaerr=0.4


x1minc=[]
x2minc=[]

x1maxc=[]
x2maxc=[]


%[c1g,c2g,c3g,c4g,c5g,c6g,c7g,d,d,d,d,d,d,d,d]=textread('datasetbig.csv', '%s %s %s %f %f %f %f %s %s %s %s %s %s %s %s' );
% read expiration dates
fid = fopen('ExpirationDates.txt','rt'); 
indata = textscan(fid, '%s %s', 'HeaderLines',1); 
fclose(fid); 
expdata = indata{1,1};
[ye, me, de] = datevec(expdata,'dd-mmm-yy');

% read DAX
[a1,a2]=textread('dax_index.txt','%s %f');
dax=horzcat(a2);
date_dax=datenum(a1,'dd.mm.yyyy'); 


% read interest rate
[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16]=textread('IRs.dat','%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
date_ir=datenum(b1,'dd.mm.yyyy');
IR=horzcat(b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

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
%[c111,c112,c113,c114,c115,c116,c117]=textread('C_2012.txt', '%s %s %s  %f %f %f %f');

c1g=[c11;c21;c31;c41;c51;c61;c71;c81;c91;c101];
c2g=[c12;c22;c32;c42;c52;c62;c72;c82;c92;c102];
c3g=[c13;c23;c33;c43;c53;c63;c73;c83;c93;c103];
c4g=[c14;c24;c34;c44;c54;c64;c74;c84;c94;c104];
c5g=[c15;c25;c35;c45;c55;c65;c75;c85;c95;c105];
c6g=[c16;c26;c36;c46;c56;c66;c76;c86;c96;c106];
c7g=[c17;c27;c37;c47;c57;c67;c77;c87;c97;c107];

%c1g=[c11];
%c2g=[c12];
%c3g=[c13];
%c4g=[c14];
%c5g=[c15];
%c6g=[c16];
%c7g=[c17];

days=unique(c2g)

N=length( days )

%N=100

%N=50
L=2
ino=normrnd(0,1,N+2,1);   
st(1)=100;  
cgridx=[]
  
cgridy=[]

vdaxc=[]
odax=[]

daysa=[]



fid = fopen('ExpirationDates.txt','rt'); 
indata = textscan(fid, '%s %s', 'HeaderLines',1); 
fclose(fid); 
expdata = indata{1,1};
[ye, me, de] = datevec(expdata,'dd-mmm-yy');

day=datenum(c2g, 'dd/mm/yyyy');

for d=2:(N+1)

d 

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

% read expiration dates


% compute maturities
clear tauc
clear aa

for i=1:length(cd)
aa(i)=find(ye==c5(i) & me==c4(i));

tauc(i)=(datenum(expdata(aa(i)),'dd-mmm-yy')-refdmy_c(i,1))/365;

end

% read and find stock price
%[a1,a2]=textread('dax_index.dat','%s %f');
%dax=horzcat(a2);
%date_dax=datenum(a1,'dd/mm/yyyy');

dday=find(date_dax==numdate);
s=dax(dday);

if(s>0)
    
daysa=[daysa, dday ];

[a1v,a2v]=textread('vdax.txt','%s %f');
vdax=horzcat(a2v);
date_vdax=datenum(a1v,'dd.mm.yyyy');

vdday=find(date_vdax==numdate);
sv=vdax(vdday);

odax=[odax,s];

vdaxc=[vdaxc, sv ];

% read and compute interest rate
[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16]=textread('IRs.dat','%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
date_ir=datenum(b1,'dd.mm.yyyy');
IR=horzcat(b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

rday=find(date_ir==numdate);
R=IR(rday,:);

T=[7/365 14/365 21/365 1/12 2/12 3/12 4/12 5/12 6/12 7/12 8/12 9/12 10/12 11/12 1]; % maturity for each IR





rc=interp1(T,R,tauc,'linear','extrap');

cdd=[cd(:,1) cd(:,4:5)/s rc'/100 tauc']; %date strike option ir mat

%cdd=[cd(:,1) cd(:,4)/s cd(:,5)/s.*exp(rc'.*tauc') rc'/100 tauc']; %date strike option ir mat


cdd(find(tauc>1),:)=[];
% plot data by maturities

%Sorting maturity
[cdd(:,5), Im] = sort( cdd(:,5) );
cdd(:,2)=cdd(Im,2);
cdd(:,3)=cdd(Im,3);

%sorting Moneyness as well
matu=unique(cdd(:,5));

%XXXX

for(l=1:length(matu) )
   Ik=find(cdd(:,5)==matu(l) );
   [cdd(Ik,2),Is] =sort(cdd(Ik,2));
   
    cdd(Ik,3)=cdd((min(Ik)-1)+Is,3);
end

%F = TriScatteredInterp(cdd(:,5),cdd(:,2),cdd(:,3), 'nearest');   %%%%Interpolate Data for FPCA
F = TriScatteredInterp(cdd(:,5),cdd(:,2),cdd(:,3), 'linear');   %%%%Interpolate Data for FPCA
%F = TriScatteredInterp(cdd(:,5),cdd(:,2),cdd(:,3));   %%%%Interpolate Data for FPCA
Fc{d-1}=F; %%Store interpolated curves

%F2 = TriScatteredInterp(cdd(:,5),cdd(:,2),cdd(:,3), 'nearest');   %%%%Interpolate Data for FPCA
%Fc2{d-1}=F2; %%Store interpolated curves



%%%find smallest common grid...
x1minc=[ x1minc min(cdd(:,5)) ];
x2minc=[ x2minc  min(cdd(:,2)) ];

x1maxc=[x1maxc max(cdd(:,5))];
x2maxc=[ x2maxc max(cdd(:,2)) ];

cgridx=union(cgridx,cdd(:,5));
cgridy=union(cgridy,cdd(:,2));
c5unil{d-1}= cdd(:,5);   %%Maturity
c2unil{d-1}=cdd(:,2);    %%Moneyness
c3unil{d-1}=cdd(:,3);    %%Call price
end
end
%%%READ REAL DATA END

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

Tcum=9999999
Tmon=9999
Tmat=9999
parfor(i=1:N)
Tmon=min( Tmon,length(unique( cell2mat( c2unil(i) ) )));
Tmat=min( Tmat, length(unique( cell2mat( c5unil(i) ) )));  
Tcum=min(Tcum, length(unique( cell2mat( c2unil(i) ) ))*length(unique( cell2mat( c5unil(i) ) )))
end

x1min   =median( x1minc );
x2min   =median( x2minc );
x1max   =median( x1maxc );
x2max   =1.4 %median( x2maxc );


%x1min   =max( x1minc );
%x2min   =max( x2minc );
%x1max   =min( x1maxc );
%x2max   =min( x2maxc );

my      =x1min + ( x1max-x1min ).*rand( 512 , 1 );
mx      =x2min + ( x2max-x2min ).*rand( 512 , 1 );



%[hX2b, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(1,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unil,N,Tmon,Tmat,mx,my,1,0,'cpu');


sigma=zeros(N,1);
parfor i=1:N
sigma(i)=varaince([cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) )],cell2mat(  c3unil(i) ) );
end

%sigma=ones(N,1)*mean(sigma);

%[hX2b, V2b, loadsb,Meansmo2bb,Db]=fpca1_gpu(3,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unil,N,mx,my,1,'cpu',sigma,1);

Xc=[]
for i=1:N
    Xc{i}= [ cell2mat(c5unil(i)),cell2mat(c2unil(i)) ]
end

%[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(3,Fc,5,[min(x1minc),min(x2minc)],[max(x1maxc),max(x2maxc)],Xc,c3unil,0,[mx,my],1,sigma);
[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(3,Fc,5,[min(x1minc),min(x2minc)],[max(x1maxc),max(x2maxc)],Xc,c3unil,0,[my,mx],1,sigma,1);


%[hX2b, V2b, loadsb,Meansmo2bb,Db,mx,my]        =fpca1app(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,0,N,Tmon,Tmat,0,0)
%[hX2b, V2b, loadsb,Meansmo2bb,Db,mx,my]        =fpca1a(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,0,N,Tmon,Tmat,0,0);
%[hX2b, V2b, loadsb,Meansmo2bb,Db,mx,my]        =fpca1app(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,0,N,Tmon,Tmat,0,0)
%[hX2bx, V2bx, loadsbx,Meansmo2bbx,Dbx,mxx,myx]        =fpca1app(2,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,0,N,Tmon,Tmat,0,0);

scatter3(my,mx ,hX2b(:,1) )   %%%Plot first density
scatter3(mx,my ,Meansmo2bb )   %%%Plot Meancurve

L=3
for i=1:3
hold on
scatter3(mx,my ,V2b(:,i) )   %%%Plot eigenfunctions
end

 plot(sigma(Is))
%%%optimal projection to vdax

[dd, Is]=sort(daysa);
vdaxm= vdaxc - mean( vdaxc);
fitlm(loadsb(1:L,Is)' ,vdaxm(Is) )

Lg=loadsb(1:L,Is)'
a=(Lg'*Lg)^(-1)*Lg'*vdaxm(Is)'      %%compute the projection

hold off
plot(Lg*a)
hold on
plot(vdaxm(Is))

%plot(Lg*a-vdaxm(Is)')

%%%Score of second component looks strange
for i=1:3
hold on
plot(Lg(:,i)*a(i) )  %%%Plot scores
end

%Maybe drop second comp. and state that it is noise. 1st and 3rd looks very
%fine though :)


%%ARMIA MODEL 

Mdl = arima('MA',{0.5,-0.3},'SMA',0.4,'SMALags',12,...
		'Constant',0.04,'Variance',0.2);
%rng(200);
%Y = simulate(Mdl,130);
farc=30
start=1900
vec=1
Lgd=diff(Lg)
ToEstMdl = arima('MALags',1:2,'SMALags',12);
EstMdl = estimate(ToEstMdl,Lgd(1:start,vec));


[YF YMSE] = forecast(EstMdl,farc,'Y0',Lgd(1:start,vec));
figure
h1 = plot(start-100:start+100,Lgd(start-100:start+100,vec),'Color',[.7,.7,.7]);
hold on
h2 = plot((start+1):(start+farc),YF,'b','LineWidth',2);
h3 = plot((start+1):(start+farc),YF + 1.96*sqrt(YMSE),'r:',...
		'LineWidth',2);
plot((start+1):(start+farc),YF - 1.96*sqrt(YMSE),'r:','LineWidth',2);
legend([h1 h2 h3],'Observed','Forecast',...
		'95% Confidence Interval','Location','NorthWest');
title(['30-Period Forecasts and Approximate 95% '...
			'Confidence Intervals'])
hold off

start=2100
vec=1
volModel = arima(1,1,1);
volFit = estimate(volModel,Lg(1:start,vec));
[Y,YMSE] = forecast(volFit,500,'Y0',Lg(1:start,vec));
lower = Y - 1.96*sqrt(YMSE);
upper = Y + 1.96*sqrt(YMSE);

figure
plot(Lg(:,vec),'Color',[.7,.7,.7]);
hold on
h1 = plot( (start+1):(start+500),lower,'r:','LineWidth',2);
plot( (start+1):(start+500),upper,'r:','LineWidth',2)
h2 = plot( (start+1):(start+500),Y,'k','LineWidth',2);
legend([h1 h2],'95% Interval','Forecast',...
	     'Location','NorthWest')
title('First component Forecast')
hold off