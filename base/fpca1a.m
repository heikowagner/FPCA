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

function [hX2r, V2a, loadsa,Meansmo2b,Da,mx,my] = fpca1(L,Fc,x1minc,x2minc,x1maxc,x2maxc,cgridx,cgridy,c5unil,c2unil,c3unil,c3unilr,N,Tmon,Tmat,mx,my,method)
%%%To compute M we need a random Grid
x1min   =max( x1minc );
x2min   =max( x2minc );
x1max   =min( x1maxc );
x2max   =min( x2maxc );
if(mx==0)
my      =x1min + ( x1max-x1min ).*rand( Tmon*Tmat , 1 );
mx      =x2min + ( x2max-x2min ).*rand( Tmon*Tmat , 1 );
end
densy   =x2max-x2min;
densx   =x1max-x1min;

%Fit the data to a common grid
Yint    =[];
for i=1:N   
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




sigmae2=[];
parfor i=1:N
  
Umon= unique(cell2mat(c5unil(i)));
dummy2=cell2mat(c3unil(i));
dummy3=cell2mat(c2unil(i));
est2=[];
we=sqrt(1.5)^(-1)*[-0.5,1,-0.5];

%sum(we.^2);

for j=1:length(Umon)
  dummy21=dummy2(Umon(j)==cell2mat(c5unil(i)));
  dummy31=dummy3(Umon(j)==cell2mat(c5unil(i)));
  [d,I]= sort(dummy31);
  est     =( we(1)*dummy21(I(1:(end-2)),: ) + we(2)*dummy21(I(2:(end-1)),: ) + we(3)*dummy21(I(3:(end)),: ));
  est2=[est2 mean(est.^2)];
end
sigmae2 =[sigmae2 mean(est2)];
end
sigmae  =sigmae2 ; 

%sigmae  =repmat(0.00000005,N,1)'

Lppd4x=[mx*0, mx*0, mx*0, 0*6*mx.^0, 24*mx.^0, 120*mx.^1, my*0, my*0, my*0, 24*my*0];
Lppd4y=[mx*0, mx*0, mx*0, 6*mx*0, 24*mx*0, 120*mx*0, my*0, my*0, 0*6*my.^0, 24*my.^0];

h1m=[];
h2m=[];
parfor i=1:N
       
%%%Estimate optimal individual smoothing parameter
%%Construct polynomial estimator for h1

xachs   =cell2mat(  c2unil(i) );
yachs   =cell2mat(  c5unil(i) );
Lp      =[xachs.^0, xachs, xachs.^2, xachs.^3, xachs.^4, xachs.^5, yachs.^1, yachs.^2, yachs.^3, yachs.^4];
Vxre    =cell2mat(  c3unil(i) );
scores  =(Lp'*Lp)^(-1)*Lp'*Vxre;
Lppm    =[mx.^0, mx, mx.^2, mx.^3, mx.^4, mx.^5, my.^1, my.^2, my.^3, my.^4];
Vd3x    =Lppd4x*scores;
Vd3y    =Lppd4y*scores;
%h1i     =(1/(Tmat))*((sigmae(i))/(densx*mean(Vd3x.^2)));
%h2i     =(1/(Tmon*Tmat))*((sigmae(i))/( densy*mean(Vd3y.^2)));
h2i=(max(cell2mat(  c2unil(i) ))-min(cell2mat(  c2unil(i) )))   *length(unique(cell2mat(  c2unil(i) ) ));
h1i=(max(cell2mat(  c5unil(i) ))-min(cell2mat(  c5unil(i) )))   *length(unique(cell2mat(  c5unil(i) ) ));
h1m     =[h1m h1i];
h2m     =[h2m h2i];
end

W=[]
Xsmo2m  =[];

parfor i=1:N
   if method==1
[XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,h1m(i)^(-1/3),h2m(i)^(-1/3),1,'Epan'  ); %%Estimate loc poly
Xsmo2m=[Xsmo2m XmiS0m] ;
W= [W W0];
   else
[XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,h1m(i)^(-1/10),h2m(i)^(-1/10),7,'Epan'  ); %%Estimate loc poly
Xsmo2m=[Xsmo2m XmiS2m] ;   
W= [W W2];
   end
   
end

% i=1
% h1i     =(1/(Tmon))^(1/10)*((sigmae(i))/(densx*mean(Vd3x.^2)))^(1/10);
% [XmiS0m XmiS1m XmiS2m XmiS3m, W0, W1, W2, W3]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,5*h1i^(10/3) ,0.3*h2m(i)^(10/3),3  ); %%Estimate loc poly
% scatter3(cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ))
% hold on
%scatter3(my,mx,XmiS0m)
%scatter3(my,mx,W0)

Xsmo2m  =real(Xsmo2m);

%Md
Xa    =Xsmo2m % -repmat(mean(Xsmo2m,2),[1 N]);
Muc   =(1/size(Xa,1))*Xa'*Xa - diag(mean(W).*sigmae) ;
M     =Muc - (1/size(Xa,1))*Xa'*repmat(mean(Xsmo2m,2),[1 N]) - (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*Xa + (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*repmat(mean(Xsmo2m,2),[1 N]);

[Vduala, Da]    =eig(M); 
Da              =Da(Da>0);
Da2             =Da;
Vduala          =Vduala(:,Da>0);
Da              =diag(Da);

h1a=[];
h2a=[];

T  =Tmon*Tmat;
Lp=[mx.^0, mx, mx.^2, mx.^3, mx.^4, mx.^5, my.^1, my.^2, my.^3, my.^4];
scores=(Lp'*Lp)^(-1)*Lp'*X;

regsx=Lppd4x*scores;
regsy=Lppd4y*scores;

for i=1:(L+1)
%h1a= [h1a (1/Tmat)^(1/10)* (mean(sigmae)/( densx*mean( (regsy*Vduala(:,i)).^2 ) )  )^(1/10)];
%h2a= [h2a (1/T)^(1/10)* (mean(sigmae)/( densy*mean( (regsx*Vduala(:,i)).^2    ) )  )^(1/10)];
h1a= [h1a (1/Tmat)^(1/10)* (mean(sigmae)/( densx*sum( (regsy*Vduala(:,i)).^2 ) )  )^(1/10)];
h2a= [h2a (1/T)^(1/10)* (mean(sigmae)/( densy*sum( (regsx*Vduala(:,i)).^2    ) )  )^(1/10)];
end

% h1a     =0.5*mean(h1a)
% h2a     =2*mean(h2a)
h1a     =mean(h1a);
h2a     =mean(h2a);

Xsmoa   =[] ;   
Xsmo2a  =[] ;
parfor i=1:N
	[XmiS0a XmiS1a XmiS2a XmiS3a]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) , cell2mat(  c3unil(i) ),my ,mx ,h1a, Cdp *h2a ,3,'Gauss' ); %%Estimate loc poly
	Xsmoa   =[Xsmoa XmiS0a];
	Xsmo2a  =[Xsmo2a XmiS2a];
end

Meansmo2b   =mean(Xsmo2a')';                %Mean curve
V2a         =(Xsmo2a )*Vduala*sqrt(Da)^-1;	%Eigenvectors 2nd derivative smoothed estimate
loadsa      = sqrt(Da)*Vduala';             %Loadings

%%%Final estimator
hX2r        =repmat(Meansmo2b,[1 N])+ V2a(:,1:L+1)*loadsa(1:L+1,:);
end
% 
% scatter3(my,mx,hX2r(:,2))
% hold on
% scatter3(my,mx,realD(:,2))

%scatter3(my,mx, Xsmo2a(:,2))
