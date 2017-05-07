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
% Input:      
%L          - Number of Dimensions
%der        - g dimensional vector of derivatives 
%method     - 0 for M(0) method and d for M(d) method 
%Xc         - g dimensional joined Grid x
%X          - g dimensional domain
%Y          - Obervation with error
%x          - Output grid (optinal)
%sigma      - Nx1 known sigma for diagonal correction (optional)
%H          - gxg Bandwidth matrix to derive gamma(optional) 
%xminc      - g dimensional vector of Min x values
%xmaxc      - g dimensional vector of Min x values


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
% Author:      Heiko Wagner, 2017/01/25
% ------------------------------------------------------------------------------ 

function [hX2r, V2a, loadsa,Meansmo2b,Eigenvals,x] = fpca1(L,der,method,Xc,Y,x,sigma,H,xminc,xmaxc, kernel, kernel2)

%%Internal grid used for integration
g= size(cell2mat(Xc(1)),2);             %Dimension of the Domain
N=size(Y,2);                            %Number of curves
Tint=size(x,1);                         %Number of points at this grid


%%Construct internal random Grid used for integration
xmin   = xminc ;
xmax   = xmaxc ;
if(Tint<2)
        'No output grid chosen, Random grid with 512 will is construced'
        Tint=512;
        x=repmat(xmin ,Tint,1) + repmat( xmax-xmin ,Tint,1).*rand( Tint , g ); 
end

%%Rough estimate of the density, used for Bandwidth selection
dens=zeros(g,1);
for i=1:g
    dens(i)   =1/( xmax(i)-xmin(i) );
end

%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE

%%%Estimate hat(sigma)

if(sigma==0)
    'No sigma chosen, sigma is estimated'
    sigma=zeros(N,1);
    parfor i=1:N
        sigma(i)=varaince(cell2mat( Xc(i) ),cell2mat(  Y(i) ) );
    end
end


%%Derive very simple bw estimator for M
hm=ones(N,g);
Tobs=zeros(1,g);
for i=1:N
    dummy=cell2mat(  Xc(i) );
        for j=1:g
            Tobs(j)=size(unique(dummy(:,j)),1);
        end
        %hm(i,:)= sigma(i)./(1/range(dummy).* Tobs );
        hm(i,:)= (1./range(dummy).* Tobs );

end


if(method=='M_0')
    'M_0 Method'
    derM=zeros(1,g);
else
    'M_d Method'
    derM=der;
end

%%Choice of parameters based on Remark 1
rho=max(sum(derM)+1, ceil(g/2-1+ 3*sum(round((sum(derM)-1)/g))));

%h_ord= 0.5*( (1/(2*(rho+1- sum(derM)) ))+ (1/(g+ sum(derM ) )));

h_ord=-1/(2*(rho+1- sum(derM)) );
%h_ord=(1/(g+ sum(derM ) ));


%h_ord= (1/(2*(rho+1- sum(derM)) ));
%h_ord=-(1/2);

%%Find position of the derivative of multiloc output
pos=find( ismember(permn(0:sum(der),g)  ,derM,'rows'),1);
pos
%Smooth curves for M Matrix
W=[];
Xsmo2m  =[];
parfor i=1:N
    [W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,diag( hm(i,:).^h_ord),rho,kernel2  ); %%Estimate loc poly           
    %[W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,diag( hm(i,:).^h_ord),rho,'Gauss'  ); %%Estimate loc poly           
    Xsmo2m=[Xsmo2m fit_M(:,pos)] ;
    W= [W, W_M(:,pos)];
end


Xa    =Xsmo2m - repmat(mean(Xsmo2m,2),[1 N]); 
%plot(x,Xa(:,1:3),'.')
M   =(prod(range(x))/size(Xa,1))*(Xa'*Xa - diag(mean(W).*sigma')) ;
%M     =( Muc - (1/size(Xa,1))*Xa'*repmat(mean(Xsmo2m,2),[1 N]) - (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*Xa + (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*repmat(mean(Xsmo2m,2),[1 N])); %%Mean correction

%M=(prod(range(x))/size(Xa,1))*Xa'*Xa
%plot(mean(Xsmo2m,2))
%plot(Xa(:,37))
plot(Xa)
%imagesc(Xa'*Xa) 

%plot(x,Xa)

%%Estimation based in Eigenvaluedecomposition of the dual matrix
%[Vduala, Da]    =eig(1/N*M); 
[Vduala, Da]    =eig(M); 
%[Da,Ind]        =sort(Da(Da>0),'descend');
Da              =diag(Da);
[Da,Ind]        =sort(Da,'descend');
Da2             =diag(Da(Da>0));
%Vduala          =Vduala(:,Da>0);
Vduala          =Vduala(:,Ind);
L= min( (size(Da2,1)-1),L ); %%drop negative eigenvalues





if(H==0)
    'No BW-Matrix, Bandwidth is estimated'
    order= floor( floor(2*(sum(der)+1))/2)
    Lpx=zeros(size(x,1),1);
    
dummy=permn(0:order,g);
k=dummy( sum(dummy,2)<order+1,: );
pos2=find( ismember(k  ,der,'rows'),1);
%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE


int=6*(1:500)/500 -3 ;
int=int';
eqkernel=equivkernel(int,order,kernel);

Cdp=zeros(1,g);
for l=1:g
Cdp(l)=factorial(order+1)^2*(2*der(l)+1)*range(int)*mean(eqkernel(der(l)+1,:).^2)/(2*(order+1-der(l) )*range(int)^2*mean(eqkernel(der(l)+1,:).*int'.^(order+1) )^(2*g) ) ;
end

    for j=1:(order+2)
        Lpx      =[Lpx, factorial(j)*x.^max(0,j-order )*min(1, max(0,j-order+1 ) )  ]; 
    end
    
    scores=[];
    parfor i=1:N
        xachs   =cell2mat(  Xc(i) );
        Lp=ones(size(xachs,1),1);
        for j=1:(order+2)
            Lp      =[Lp, xachs.^j];
        end
        Vxre    =cell2mat(  Y(i) );
        scores  =[scores (Lp'*Lp)^(-1)*Lp'*Vxre];
    end
    regsx=Lpx*scores;
    ha=zeros(g,(L));
    for i=1:L

       % ha(:,i)= mean(sigma)./( Tobs'.*( dens.*mean( (regsx*Vduala(:,i)).^2  )));
       ha(:,i)= sigma'*(Vduala(:,i).^2)./( mean(Tobs)*( dens.*range(x)'.*mean( ( regsx*Vduala(:,i) ).^2  )));
      % ha(:,i)= sigma'*(Vduala(:,i).^2)./( mean(Tobs)*( dens.*N));
      
    end
    H=g*diag( (mean(ha')'.*Cdp').^(1/(2*order+2+g) ) );   
end

Xsmoa   =[] ;   
Xsmo2a  =[] ;
parfor i=1:N
   % [W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,H,order,'Gauss'  );
     [W_dm,fit]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,H,order,kernel  );
    Xsmo2a  =[Xsmo2a fit(:,pos2)];
end

Meansmo2b   =mean(Xsmo2a')';                %Mean curve
V2a         =(Xsmo2a-repmat(Meansmo2b,[1 N]) ) *Vduala*diag((Da).^(-1/2));	%Eigenvectors 2nd derivative smoothed estimate
loadsa      = diag((Da).^(1/2))*Vduala';             %Loadings
Eigenvals=Da/N;

%%%Final estimator
hX2r        =repmat(Meansmo2b,[1 N])+ V2a(:,1:L+1)*loadsa(1:L+1,:);
%hX2r        =repmat(Meansmo2b,[1 N])+ (Xsmo2a-repmat(Meansmo2b,[1 N]) ) *Vduala(:,1:L)*Vduala(:,1:L)';
%hX2r        =repmat(Meansmo2b,[1 N])+ Xa *Vduala(:,1:L)*Vduala(:,1:L)';

end