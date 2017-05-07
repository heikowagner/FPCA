
function [Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,realD,mx,my]=simulate(N,ta,r,mu,sigma,sigmaerr,k,Tmon,Tmat)   


%Simulate data
clear TAU ino st  p X C Q Y h1 h2
TAU  =(2:6)/7
randn('state',k+6)
ino  =randn(1,N+2)'
st(1)=100;  

% simulation  of stock price
for i=2:(N+1)
    st(i)=st(i-1)*exp((mu-0.5*sigma^2)*ta+sigma*sqrt(ta)*ino(i));  
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

m   =[0.4 0.7 0.1];
s   =[0.5 0.3 0.3]; % standard deviation
rand('state',k+1)
p=abs( randn(3,N+1 )  );
p=p./repmat(sum(p),3,1);
for i=2:(N+1)
TAU=((1:Tmat)+ rand(Tmat,1)')/(Tmat+1);
    Xj  =[];
    Yj  =[];
    Cj  =[];
    qj  =[];
    tauj=[];
    
    for j=1:Tmat  % number of maturities per day
        tau     =TAU(j)+i*12-floor(TAU(j)+i*12);
        x       =st(i)*sort( rand(Tmon,1))*2; %Random Grid
    
        for l=1:3
            qu(l,:)=1./(x*sqrt(2*pi*s(l)^2*tau)).*exp(-1/2*((log(x/st(i))- (m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau))).^2);    
            Ca(l,:)=exp(-r*tau)*st(i)*exp(m(l)*tau)*normcdf((log(st(i)./x)+(m(l)+s(l)^2/2)*tau)/(s(l)*sqrt(tau)))-exp(-r*tau)*x.*normcdf((log(st(i)./x)+(m(l)-s(l)^2/2)*tau)/(s(l)*sqrt(tau)));
        end

    Ca(isnan(Ca)) = 0 ;
    X(j,:,i)=x;   %%Observation Grid K
    Q(j,:,i)=p(:,i)'*qu;   %%Value of risk neutral density
    C(j,:,i)=p(:,i)'*Ca;
    Y(j,:,i)=C(j,:,i)+sigmaerr*randn(1,Tmon );    %%Obervation with error
    Xj=[Xj ; X(j,:,i)'];  % strike
    Yj=[Yj ; Y(j,:,i)'];  % calls
    Cj=[Cj ; C(j,:,i)'];  % calls wo error
    tauj=[tauj ; tau*ones(size(C,2),1)]; % moneyness
    qj=[qj ; Q(j,:,i)'];
    end   
    
cdd=[i*ones(length(Xj),1) Xj/st(i) Yj/st(i).*exp(r*tauj) r*ones(length(Xj),1) tauj ]; %date strike option ir mat

%%%%%Simulation finish
F           =TriScatteredInterp(cdd(:,2),cdd(:,5),cdd(:,3), 'nearest');  
Fc{i-1}     =F;                                                       %% Store interpolated curves

%%%find smallest common grid...
x1minc      =[x1minc min(cdd(:,5))];        % mat
x2minc      =[x2minc  min(cdd(:,2))];       % mon
x1maxc      =[x1maxc max(cdd(:,5))];
x2maxc      =[x2maxc max(cdd(:,2))];
cgridx      =union(cgridx,cdd(:,5));
cgridy      =union(cgridy,cdd(:,2));
c5unil{i-1} =cdd(:,5);                      % maturity
c2unil{i-1} =cdd(:,2);                      % strike/s*exp(rtau)
c3unil{i-1} =cdd(:,3);                      % call/s
c3unilr{i-1}=Cj/st(i).*exp(r*tauj); 
end
%%%%%GENERATING SIMULATION ENDE

%%%true curves
X2=X

%x1min   =max( x1minc );
%x2min   =max( x2minc );
%x1max   =min( x1maxc );
%x2max   =min( x2maxc );

x1min   =min( x1minc );
x2min   =min( x2minc );
x1max   =max( x1maxc );
x2max   =max( x2maxc );

my      =x1min + ( x1max-x1min ).*rand( Tmon*Tmat , 1 );
mx      =x2min + ( x2max-x2min ).*rand( Tmon*Tmat , 1 );
X=mx;
T=my;
for l=1:3
QQ(l,:)=1./(X.*sqrt(2*pi*s(l)^2.*T)).*exp(-1/2*((log(X)- (m(l)-s(l)^2/2)*T)./(s(l)*sqrt(T))).^2);    
end
realD=QQ'*p(:,2:end);

[Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr]
end
