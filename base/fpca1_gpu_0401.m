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
%xminc     - g dimensional vector of Min x values
%xmaxc     - g dimensional vector of Min x values
%cgridx     - g dimensional joined Grid x
%X          - g dimensional domain
%Y           - Obervationw with error
%Ytrue       - Obervations wo. error
%N          - Number of Days  
%x          -Output grid
%method     - 0 for M(0) method and d for M(d) method
%comp       -Use gpu or cpu
%sigma       - known sigma for diagonal correctin
%app

  
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
% Author:      Heiko Wagner, 2016/12/22
% ------------------------------------------------------------------------------ 

function [hX2r, V2a, loadsa,Meansmo2b,Da,x] = fpca1(L,Fc,d,xminc,xmaxc,cgridx,Xc,Y,Ytrue,x,method,comp,sigma)

% 
% L=2
% Fc=Fc
% xminc=[x1minc',x2minc']
% xmaxc=[x1maxc',x2maxc']
% cgrid=[cgridx,cgridy]
% Xc=Xc
% Y=c3unil
% Ytrue=c3unilr
% x=[mx,my]
% method=1
% comp='cpu'
% app=0
% sigma=0
Tint=512;

g=size(cell2mat(Xc(1)),2)
N=size(Y,2)
L=L-1;
mxe=x;
%%%To compute M we need a random Grid
%xmin   =max( xminc );
%xmax   =min( xmaxc );

xmin   = xminc ;
xmax   = xmaxc ;

 x      =repmat(xmin ,Tint,1) + repmat( xmax-xmin ,Tint,1).*rand( Tint , g );
 dens=zeros(g,1)
for i=1:g
    dens(i)   =xmax(i)-xmin(i);
end

%Fit the data to a common 
T    =zeros(N,g);

Yint=[];
parfor i=1:N  
    for m=1:g
    T(i,m)=length(unique(cell2mat(  Xc(m) ) ));
    end
    
    ddum=Fc{i}( x );
    Yint=[Yint ddum(:)];
end
X               =Yint- repmat( mean( Yint , 2 ), [1 N]);
X( isnan(X) )   =0;

%%%Grid construction END
 
%%%Estimate h0 in each direction

%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE
Cdp=1  %1.0006;   %%Constant for derivatives taken from example from Fan 1996
%%%Estimate hat(sigma)

if(sigma==0 )
sigma=zeros(N,1);
parfor i=1:N
    sigma(i)=varaince(cell2mat( Xc(i) ),cell2mat(  Y(i) ) );
end
end
sigmae  =sigma;


%%Derive very simple bw estimator for M
hm=ones(N,g);
parfor i=1:N
    hm(i,:)=  (max(cell2mat(  Xc(i) ))-(min(cell2mat(  Xc(i) )))).* size(cell2mat(  Xc(i) ),1)  ;
end


W=[];
Xsmo2m  =[];

for i=1:N
    
   if(method==1)
       if comp=='gpu'
            [W_M,fit_M]=multiloc( gpuArray(cell2mat( Xc(i) )) , gpuArray(cell2mat(  Y(i) )),gpuArray(x) ,gpuArray(hm(i,:).^(-1/3)),1,'Epan'  ); %%Estimate loc poly 
           else 
            [W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,hm(i,:).^(-1/3),1,'Epan'  ); %%Estimate loc poly           
       end
Xsmo2m=[Xsmo2m fit_M(:,1)] ;
W= [W, W_M(:,1)];
   else
              if comp=='gpu' 
                   [W_M,fit_M]=multiloc( gpuArray(cell2mat( Xc(i) )) , gpuArray(cell2mat(  Y(i) )),gpuArray(x) ,gpuArray(hm(i,:).^(-1/10)),7,'Epan'  );
              else
                
                   [W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,hm(i,:).^(-1/10),7,'Epan'  );
              end              
Xsmo2m=[Xsmo2m fit_M(:,d)] ;   
W= [W, W_M(:,d)];
   end
 
end

if comp=='gpu'
Xsmo2m=gpuArray(Xsmo2m);
end
size(sigmae)
size(mean(W))
size(diag(mean(W).*sigmae'))

Xa    =Xsmo2m; % -repmat(mean(Xsmo2m,2),[1 N]);
size(Xa'*Xa)
Muc   =(1/size(Xa,1))*(Xa'*Xa - diag(mean(W).*sigmae')) ;
M     =Muc - (1/size(Xa,1))*Xa'*repmat(mean(Xsmo2m,2),[1 N]) - (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*Xa + (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*repmat(mean(Xsmo2m,2),[1 N]);


[Vduala, Da]    =eig(M); 
[Da,Ind]        =sort(Da(Da>0),'descend');
Da2             =Da
Vduala          =Vduala(:,Da>0);
Vduala          =Vduala(:,Ind);
Da              =diag(Da);

L= min( (size(Da2,1)-1),L );

Lpx      =[ones(size(x,1),1).^0, x.^0, x.^0, x.^0, 4*3*2*1*x]; %%has to be replaced with order+ 1 derivative of xachs
scores=[];
for i=1:N
xachs   =cell2mat(  Xc(i) );
Lp      =[ones(size(xachs,1),1), xachs, xachs.^2, xachs.^3, xachs.^4];
Vxre    =cell2mat(  Y(i) );
scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
end

regsx=Lpx*scores;

%regsx1=Lp*scores;

%scatter3(x(:,1), x(:,2), regsx(:,i) )
%hold on
% scatter3(xachs(:,1), xachs(:,2), cell2mat( Y(i))  )
% hold on 
% scatter3(xachs(:,1), xachs(:,2), regsx1(:,i)  )


ha=zeros(2,(L+1));
Tobs=size(xachs,1)
for i=1:(L+1)
ha(:,i)= mean(sigmae)./( Tobs*( dens.*mean( (regsx*Vduala(:,i)).^2  )));

end

h1a=mean(ha)
% h1a     =0.5*h1a
% h2a     =2*mean(h2a)

%h1a     =repmat( mean(ha),N,1);

Xsmoa   =[] ;   
Xsmo2a  =[] ;
parfor i=1:N
    if comp=='gpu'
        %[XmiS0a XmiS1a XmiS2a XmiS3a]=multiloc( gpuArray(cell2mat(  c5unil(i) )),gpuArray(cell2mat(  c2unil(i) )) , gpuArray(cell2mat(  c3unil(i) )),my ,mx ,h1a, Cdp *h2a ,3,'Epan' ); %%Estimate loc poly
        [W_dm,fit]=multiloc( gpuArray(cell2mat(  c5unil(i) )),gpuArray(cell2mat(  c2unil(i) )) , gpuArray(cell2mat(  c3unil(i) )),mye ,mxe ,[h1a(1)^(1/10), Cdp *h1a(2)^(1/10)] ,3,'Gauss' ); %%Estimate loc poly
    
    else
        %[XmiS0a XmiS1a XmiS2a XmiS3a]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) ,cell2mat(  c3unil(i) ),my ,mx ,h1a, Cdp *h2a ,3,'Epan' ); %%Estimate loc poly
          [W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,[h1a(1)^(1/10), Cdp *h1a(2)^(1/10)],3,'Gauss'  );
        % [XmiS0a XmiS1a XmiS2a XmiS3a]=multiloc( cell2mat(  c5unil(i) ),cell2mat(  c2unil(i) ) ,cell2mat(  c3unil(i) ),mye ,mxe , [h1a(1)^(-1/10), Cdp *h2a(2)^(-1/10)] ,3,'Gauss' ); %%Estimate loc poly
    end
    %Xsmoa   =[Xsmoa XmiS0a'];
	Xsmo2a  =[Xsmo2a fit(:,d)];
end

if comp=='gpa'
%Xsmoa=gpuArray(Xsmoa);
Xsmo2a=gpuArray(Xsmo2a);
end

Meansmo2b   =mean(Xsmo2a')';                %Mean curve
V2a         =(Xsmo2a )*Vduala*sqrt(Da)^-1;	%Eigenvectors 2nd derivative smoothed estimate
loadsa      = sqrt(Da)*Vduala';             %Loadings


%%%Final estimator
hX2r        =repmat(Meansmo2b,[1 N])+ V2a(:,1:L+1)*loadsa(1:L+1,:);
end