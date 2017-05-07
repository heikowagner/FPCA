%%One dimensional simulation

function [x1min,x1max,X_out,Y,Y_true,Ytrue,Y_Dtrue,mx,YD_true]=simu1(N,ta,r,mu,sigmaerr,k,T,Tout)   

%Simulate data
clear TAU ino st  p X C Q Y h1 h2;
randn('state',k+6);
x1minc  =[];
x1maxc  =[];

%mu_l and sigma_l for mixture
m   =[0.4 0.7 0.1];
s   =[0.5 0.3 0.3]; % standard deviation


rand('state',k+1);

%%mixture weights p
p=abs( randn(3,N )  );
p=p./repmat(sum(p),3,1);

Y=zeros(T,N);
Y_true=zeros(T,N);
YD_true=zeros(T,N);
X_out=zeros(T,N);
for i=1:N
    %x=0.5+rand( 1 , T )*1.3; %Random Grid
    x=0.3+rand( 1 , T )*1.8; %Random Grid
    X_out(:,i)=x;
    tau=0.5;
    Ca=zeros(3,T);
    qu=zeros(3,T);
    for l=1:3
            qu(l,:)=1./(x.*sqrt(2*pi*s(l)^2*tau)).*exp(-1/2*((log(x)- (m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau))).^2);    
            Ca(l,:)=exp(m(l)*tau)*normcdf((log(1./x)+(m(l)+s(l)^2/2)*tau)/(s(l)*sqrt(tau)))-x.*normcdf((log(1./x)+(m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau)));
    end
    Y(:,i)=Ca'*p(:,i)+sigmaerr*randn(T,1 );
    Ytrue(:,i)=Ca'*p(:,i);
  %  YD_true(:,i)=qu'*p(:,i);
end
%%%%%GENERATING SIMULATION ENDE
%mx      =rand( Tout , 1 )*3;        
%mx      =0.5+ rand( Tout , 1 )*1.3;        
%mx      =0.5+ (1:Tout)/Tout*1.3;        
mx      =0.3+ (1:Tout)/Tout*1.8;        
X=mx;
T=tau;
for l=1:3
QQ(l,:)=1./(X.*sqrt(2*pi*s(l)^2.*T)).*exp(-1/2*((log(X)- (m(l)-s(l)^2/2)*T)./(s(l)*sqrt(T))).^2);    
end

for l=1:3
QQo(l,:)=           exp(m(l)*tau)*normcdf((log(1./X)+(m(l)+s(l)^2/2)*tau)/(s(l)*sqrt(tau)))-X.*normcdf((log(1./X)+(m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau)));
end

Y_Dtrue=QQ'*p;
Y_true=QQo'*p;
x1min=min(mx)
x1max=max(mx)
end
