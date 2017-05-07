
function [Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simulate2(N,ta,r,mu,sigma,sigmaerr,k,T,Tout)   


%Simulate data
clear TAU ino st  p X C Q Y h1 h2;
TAU  =(2:6)/7;
randn('state',k+6);
ino  =randn(1,N+2)';
st(1)=100;  

% simulation  of stock price
for i=2:(N+1)
    %st(i)=st(i-1)*exp((mu-0.5*sigma^2)*ta+sigma*sqrt(ta)*ino(i));  
    st(i)=100;
end

x1minc  =[];
x2minc  =[];
x1maxc  =[];
x2maxc  =[];
cgridx  =[];  
cgridy  =[];
Qs      =[];
Xs      =[];
taugrid =[];

%m   =[0.4 0.7 0.1];
%s   =[0.5 0.3 0.3]; % standard deviation
s   =[0.5 0.3 0.4]; % standard deviation
m   =[0.4 0.8 0.1];
    
rand('state',k+1);
p=abs( randn(3,N+1 )  );
p=p./repmat(sum(p),3,1);
for i=2:(N+1)

    Xj  =[];
    Yj  =[];
    Cj  =[];
    qj  =[];
    tauj=[];
   % plot(st)
      % number of maturities per day
        TAU= 0.2+rand( 1 , T )*0.5 ;
        tau= TAU;
        %tau     =TAU+i*12-floor(TAU+i*12);
        x       =0.5+rand( 1 , T )*1.3; %Random Grid
        %x       =sort( rand(1,1))*2; %Random Grid
      for l=1:3
            qu(l,:)=1./(x.*sqrt(2*pi*s(l)^2*tau)).*exp(-1/2*((log(x)- (m(l)-s(l)^2/2)*tau)./(s(l)*sqrt(tau))).^2);    
           % Ca(l,:)=exp(-r*tau)*st(i)*exp(m(l)*tau)*normcdf((log(st(i)./x)+(m(l)+s(l)^2/2)*tau)/(s(l)*sqrt(tau)))-exp(-r*tau)*x.*normcdf((log(st(i)./x)+(m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau)));
            Ca(l,:)=exp(m(l)*tau).*normcdf((log(1./x)+(m(l)+s(l)^2/2).*tau)./(s(l)*sqrt(tau)))-x.*normcdf((log(1./x)+(m(l)-s(l)^2/2).*tau)./(s(l)*sqrt(tau)));
            hold on
            scatter3(tau,x,qu(l,:))
        end
        
    Ca(isnan(Ca)) = 0 ;
    X(:,:,i)=x;   %%Observation Grid K
    Q(:,:,i)=p(:,i)'*qu;   %%Value of risk neutral density
    C(:,:,i)=p(:,i)'*Ca;
    Y(:,:,i)=C(:,:,i)+sigmaerr*randn(1,1 )/10;    %%Obervation with error
    Xj=[Xj ; X(:,:,i)'];  % strike
    Yj=[Yj ; Y(:,:,i)'];  % calls
    Cj=[Cj ; C(:,:,i)'];  % calls wo error
    tauj=[tauj ; tau']; % moneyness
    qj=[qj ; Q(:,:,i)'];
end   

    
%cdd=[i*ones(length(Xj),1) Xj Yj r*ones(length(Xj),1) tauj ]; %date strike option ir mat

%%%%%Simulation finish
F           =TriScatteredInterp(x',tau',Yj, 'nearest');  
Fc{i-1}     =F;                                                       %% Store interpolated curves

%%%find smallest common grid...
x1minc      =[x1minc min(tau)];        % mat
x2minc      =[x2minc  min(x)];       % mon
x1maxc      =[x1maxc max(tau)];
x2maxc      =[x2maxc max(x)];
cgridx      =union(cgridx,tau);
cgridy      =union(cgridy,x);
c5unil{i-1} =tau;                      % maturity
c2unil{i-1} =x;                      % strike/s*exp(rtau)
c3unil{i-1} =Yj;                      % call/s
c3unilr{i-1}=Cj/st(i).*exp(r*tauj); 
%%%%%GENERATING SIMULATION ENDE

%%%true curves
X2=X;

%x1min   =max( x1minc );
%x2min   =max( x2minc );
%x1max   =min( x1maxc );
%x2max   =min( x2maxc );

x1min   =min( x1minc );
x2min   =min( x2minc );
x1max   =max( x1maxc );
x2max   =max( x2maxc );

my      =0.2+rand( Tout , 1 )*0.5;
mx      =0.5+rand( Tout , 1 )*1.3;

 %TAU= rand(1,1);
 %x       =sort( rand(1,1))*2; %Random Grid
       
        
        
X=mx;
T=my;
for l=1:3
QQ(l,:)=1./(X.*sqrt(2*pi*s(l)^2.*T)).*exp(-1/2*((log(X)- (m(l)-s(l)^2/2)*T)./(s(l)*sqrt(T))).^2);    
%QQ(l,:)=1./(X.*sqrt(2*pi*s(l)^2.*T)).*exp(-1/2*((log(X/st(1))- (m(l)-s(l)^2/2)*T)./(s(l)*sqrt(T))).^2);    
%qu(l,:)=1./(x*sqrt(2*pi*s(l)^2*tau)).*exp(-1/2*((log(x/st(i))- (m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau))).^2);           
end
realD=QQ'*p(:,2:end);

[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr];
end
