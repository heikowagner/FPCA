% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: Low dimensional decomposition of derivatives for g dimensioal surfaces
%              via dual method and  automatic bandwith selection
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
%L          - Number of Dimensions
%Fc         - Triscattered intepolated curves
%dg          - Derivative to be estimated, examples: [0 0] -> 1, [1 0] ->2, [0 1] ->3, [2 0] ->4, [0 2] ->5, [2 1] ->6 etc.
%xminc      - g dimensional vector of Min x values
%xmaxc      - g dimensional vector of Min x values
%Xc         - g dimensional joined Grid x
%X          - g dimensional domain
%Y          - Obervationw with error
%Ytrue      - Obervations wo. error (optional for simulations) 
%x          - Output grid
%method     - 0 for M(0) method and d for M(d) method
%sigma      - known sigma for diagonal correction

  
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
% Author:      Heiko Wagner, 2017/01/04
% ------------------------------------------------------------------------------ 

function [hX2r, V2a, loadsa,Meansmo2b,Da,x] = fpca1(L,Fc,dg,xminc,xmaxc,Xc,Y,Ytrue,x,method,sigma,app)

%%Internal grid used for integration
Tint=size(x,1);
g= size(cell2mat(Xc(1)),2);
N=size(Y,2);
L=L-1;
mxe=x;


dummy=permn(0:3,g);             %%%Get all possible combinations
ddg=dummy(dg);  %%%Get the derivative at the dg position


%%Construct internal random Grid used for integration
xmin   = xminc ;
xmax   = xmaxc ;
if(Tint<2)
        Tint=512;
        x=repmat(xmin ,Tint,1) + repmat( xmax-xmin ,Tint,1).*rand( Tint , g ); 
end
dens=zeros(g,1);

for i=1:g
    dens(i)   =xmax(i)-xmin(i);
end

% %Fit the data to a common 
% T    =zeros(N,g);
% Yint=[];
% parfor i=1:N  
%     for m=1:g
%         T(i,m)=length(unique(cell2mat(  Xc(m) ) ));
%     end
%     ddum=Fc{i}( x );
%     Yint=[Yint ddum(:)];
% end
% X               =Yint- repmat( mean( Yint , 2 ), [1 N]);
% X( isnan(X) )   =0;

%%%Grid construction END
 
%%%Estimate h0 in each direction

%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE
Cdp=1 ; %1.0006;   %%Constant for derivatives taken from example from Fan 1996
%%%Estimate hat(sigma)

if(sigma==0)
sigma=zeros(N,1);
parfor i=1:N
    sigma(i)=varaince(cell2mat( Xc(i) ),cell2mat(  Y(i) ) );
end
end
sigmae  =sigma;


%%Derive very simple bw estimator for M
hm=ones(N,g);
Tobs=zeros(1,g);
for i=1:N
    dummy=cell2mat(  Xc(i) );
        for j=1:g
            Tobs(j)=size(unique(dummy(:,j)),1);
        end
    hm(i,:)=  (max(dummy)-min(dummy)).* Tobs;
end
mean(hm).^(-1/3)
if(app==1)
hm= repmat( mean(hm),[N 1]) ;
end


%%Choice of parameters based on Remark 1
if(method==1)
    rho=min(1, ceil(g/2-1) );
    h_ord=-(1/(2*rho + 2) );
    pos=1;
else
    rho=max(1, ceil(g/2-1+ 3*sum(round((sum(ddg)-1)/g))));
    h_ord=-(1/(2*rho + 2 - 2*sum(round((sum(ddg)-1)/g))) );
    pos=dg;
end

W=[];
Xsmo2m  =[];
for i=1:N
    [W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,hm(i,:).^h_ord,rho,'Epan'  ); %%Estimate loc poly           
    Xsmo2m=[Xsmo2m fit_M(:,pos)] ;
    W= [W, W_M(:,pos)];
end

Xa    =Xsmo2m; 
Muc   =(1/size(Xa,1))*(Xa'*Xa - diag(mean(W).*sigmae')) ;
M     =Muc - (1/size(Xa,1))*Xa'*repmat(mean(Xsmo2m,2),[1 N]) - (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*Xa + (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*repmat(mean(Xsmo2m,2),[1 N]); %%Mean correction

 
%%Estimation based in Eigenvaluedecomposition of the dual matrix
[Vduala, Da]    =eig(1/N*M); 
[Da,Ind]        =sort(Da(Da>0),'descend');
Da2             =Da;
Vduala          =Vduala(:,Da>0);
Vduala          =Vduala(:,Ind);
Da              =diag(Da);

L= min( (size(Da2,1)-1),L ); %%drop negative eigenvalues

Lpx      =[ones(size(x,1),1)*0, x*0, x*0, x*0, 24*x]; %%has to be replaced with order+ 1 derivative of xachs
scores=[];
for i=1:N
xachs   =cell2mat(  Xc(i) );
Lp      =[ones(size(xachs,1),1), xachs, xachs.^2, xachs.^3, xachs.^4];
Vxre    =cell2mat(  Y(i) );
scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
end
regsx=Lpx*scores;

ha=zeros(g,(L+1));
%Tobs=size(xachs,1)

for i=1:(L+1)
ha(:,i)= mean(sigmae)./( Tobs'.*( dens.*N*mean( (regsx*Vduala(:,i)).^2  )));
end
h1a=mean(ha');
 (regsx*Vduala(:,i))

Xsmoa   =[] ;   
Xsmo2a  =[] ;
parfor i=1:N
    [W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,h1a.^(1/(2*3+2+g) ),3,'Gauss'  );
    Xsmo2a  =[Xsmo2a fit(:,dg)];
end

Meansmo2b   =mean(Xsmo2a')';                %Mean curve
V2a         =(Xsmo2a )*Vduala*sqrt(Da)^-1;	%Eigenvectors 2nd derivative smoothed estimate
loadsa      = sqrt(Da)*Vduala';             %Loadings


%%%Final estimator
hX2r        =repmat(Meansmo2b,[1 N])+ V2a(:,1:L+1)*loadsa(1:L+1,:);
end