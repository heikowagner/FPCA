% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: Low dimensional decomposition of derivatives for surfaces
%              via method M(d) with automatic bandwith selection
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
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
% Keywords:    FPCA, Surfaces, Derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2015/08/13 
% ------------------------------------------------------------------------------ 

function [hX2r,mx,my] = individual(Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,mx,my,method,comp,sigma,order,kernel)
mxe=mx;
mye=my;
%%%To compute M we need a random Grid
x1min   =max( x1minc );
x2min   =max( x2minc );
x1max   =min( x1maxc );
x2max   =min( x2maxc );

Tint=length(cell2mat(c5unil(1)) );

%if(mx==0)
my      =x1min + ( x1max-x1min ).*rand( Tint , 1 );
mx      =x2min + ( x2max-x2min ).*rand( Tint , 1 );
%end
densy   =1/( x2max-x2min );
densx   =1/ (x1max-x1min );

%Fit the data to a common 
Yint    =[];
Tmat    =zeros(N,1);
Tmon    =zeros(N,1);

scores=[];
parfor i=1:N  
    
    xachs   =cell2mat(  c2unil(i) );
    yachs   =cell2mat(  c5unil(i) );
    Lp      =[xachs.^0, xachs, xachs.^2, xachs.^3, xachs.^4, xachs.^5, yachs.^1, yachs.^2, yachs.^3, yachs.^4];
    Vxre    =cell2mat(  c3unil(i) );
    scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
    

    Tmat(i)=length(unique(cell2mat(  c2unil(i) ) ));
    Tmon(i)=length(unique(cell2mat(  c5unil(i) ) ));
  
    ddum=Fc{i}( mx(:) , my(:) );
    Yint=[Yint ddum(:)];
end
X               =Yint- repmat( mean( Yint , 2 ), [1 N]);
X( isnan(X) )   =0;

%%%Grid construction END
 
%%%Estimate h0 in each direction

%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE
Cdp=1.0006;   %%Constant for derivatives taken from Fan 1996
%%%Estimate hat(sigma)

if(sigma==0 )
sigma=zeros(N,1);

parfor i=1:N
sigma(i)=varaince(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ),cell2mat(  c3unil(i) ) );
end
end
sigmae  =sigma;


g=2;
der=[2 0];
x=[cgridx cgridy];
Lpx=zeros(size(x,1),1);
int=6*(1:500)/500 -3 ;
int=int';
eqkernel=equivkernel(int,order,kernel);

dens=[densx, densy];


Cdp=zeros(1,g);
for l=1:g
Cdp(l)=factorial(order+1)^2*(2*der(l)+1)*range(int)*mean(eqkernel(der(l)+1,:).^2)/(2*(order+1-der(l) )*range(int)^2*mean(eqkernel(der(l)+1,:).*int'.^(order+1) )^(2*g) ) ;
end
 

    for j=1:(order+2)
        Lpx      =[Lpx, factorial(j)*x.^max(0,j-order )*min(1, max(0,j-order+1 ) )  ]; 
    end

    
    scores=[];
    for i=1:N
        xachs   =[ cell2mat(c2unil(i)) cell2mat(c5unil(i)) ];
        Lp=ones(size(xachs,1),1);
        for j=1:(order+2)
            Lp      =[Lp, xachs.^j];
        end
        Vxre    =cell2mat(  c3unil(i) );
        scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
    end

    
    regsx=Lpx*scores;
    
    
    ha=mean(sigma)./(Tint*( dens'.*range(x)'.*mean( mean( ( regsx ).^2  ))));
    H=g*diag((  ha.*Cdp').^(1/(2*order+2+g) ) )   

% % % 
% % % 
% % % h1a=[];
% % % h2a=[];
% % % 
% % % Lppd4x=[mx*0, mx*0, mx*0, 0*6*mx.^0, 24*mx.^0, 120*mx.^1, my*0, my*0, my*0, 24*my*0];
% % % Lppd4y=[mx*0, mx*0, mx*0, 6*mx*0, 24*mx*0, 120*mx*0, my*0, my*0, 0*6*my.^0, 24*my.^0];
% % % 
% % % regsx=Lppd4x*scores;
% % % regsy=Lppd4y*scores;
% % % 
% % % parfor i=1:N
% % % h1a= [h1a (1/Tmat(i))^(1/10)* (sigmae(i)/( densx* (mean(regsy(:,i) ).^2 ) )  )^(1/10)];
% % % h2a= [h2a (1/Tmon(i))^(1/10)* (sigmae(i)/( densy* (mean(regsx(:,i) ).^2    ) )  )^(1/10)];
% % % end

% h1a     =0.5*mean(h1a)
% h2a     =2*mean(h2a)
%h1a     =mean(h1a);
%h2a     =mean(h2a);

Xsmoa   =[] ;   
Xsmo2a  =[] ;
parfor i=1:N
    
        [W_M,fit]=multiloc( [cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) )] ,cell2mat(  c3unil(i) ),[mxe ,mye] ,H ,order,kernel ); %%Estimate loc poly
       % [W_M,fit_M]                  =multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,hm(i,:).^h_ord,rho,'Epan'  ); %%Estimate loc poly           
    
    Xsmoa   =[Xsmoa fit(:,1)];
	Xsmo2a  =[Xsmo2a fit(:,8)];
end


%%%Final estimator
hX2r        =Xsmo2a;
end
% 
% scatter3(my,mx,hX2r(:,2))
% hold on
% scatter3(my,mx,realD(:,2))

%scatter3(my,mx, Xsmo2a(:,2))
