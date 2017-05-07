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

function [hX2r, V2a, loadsa,Meansmo2b,Da,x] = fpca1(L,der,method,Xc,Y,x,sigma,H,xminc,xmaxc,kernel)

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
    dens(i)   =xmax(i)-xmin(i);
end

%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE
Cdp=1 ; %1.0006;   %%Constant for derivatives taken from example from Fan 1996
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
    hm(i,:)= sigma(i)/((max(dummy)-min(dummy)).* Tobs);
    %hm(i,:)= 1/((max(dummy)-min(dummy)).* Tobs);
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
h_ord= 0.5*( (1/(2*(rho+1- sum(derM)) ))+ (1/(g+ sum(derM ) )));
%h_ord= (1/(2*(rho+1- sum(derM)) ));
%h_ord=-(1/2);

%%Find position of the derivative of multiloc output
pos=find( ismember(permn(0:sum(der),g)  ,derM,'rows'),1);
pos
%Smooth curves for M Matrix
W=[];
Xsmo2m  =[];
parfor i=1:N
    [W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,diag( hm(i,:).^h_ord),rho,'Epan'  ); %%Estimate loc poly           
    %[W_M,fit_M]=multiloc( cell2mat( Xc(i) ) , cell2mat(  Y(i) ),x ,diag( hm(i,:).^h_ord),rho,'Gauss'  ); %%Estimate loc poly           
    Xsmo2m=[Xsmo2m fit_M(:,pos)] ;
    W= [W, W_M(:,pos)];
end

Xa    =Xsmo2m; 
Muc   =(1/size(Xa,1))*(Xa'*Xa - diag(mean(W).*sigma')) ;
M     =Muc - (1/size(Xa,1))*Xa'*repmat(mean(Xsmo2m,2),[1 N]) - (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*Xa + (1/size(Xa,1))*repmat(mean(Xsmo2m,2),[1 N])'*repmat(mean(Xsmo2m,2),[1 N]); %%Mean correction

plot(x,Xa)

%%Estimation based in Eigenvaluedecomposition of the dual matrix
%[Vduala, Da]    =eig(1/N*M); 
[Vduala, Da]    =eig(M); 
[Da,Ind]        =sort(Da(Da>0),'descend');
Da2             =Da;
Da2(1:2)
Vduala          =Vduala(:,Da>0);
Vduala          =Vduala(:,Ind);
Da              =diag(Da);
L= min( (size(Da2,1)-1),L ); %%drop negative eigenvalues





if(H==0)
    'No BW-Matrix, Bandwidth is estimated'
    order= (sum(der)+1);
    Lpx=zeros(size(x,1),1);
    
pos2=find( ismember(permn(0:sum(der),g)  ,der,'rows'),1);

%%%ADD DERIVAITVE CONSTANT FOR 1ST DERIVATIVE HERE

%order=2
%der=1
%pos2=der+1
int=((-300:300)/100)';
eqkernel=equivkernel(int,order,kernel);
%plot(eqkernel')
%range(int)*mean(eqkernel')

Cdp=4*factorial(order+1)^2*(2*sum(der)+1)*range(int)*mean(eqkernel(pos2,:).^2)/(2*(order+1-sum(der) )*range(int)^2*mean(int'.^(order+1).*eqkernel(pos2,:) )^2 ) ;
Cdp^(1/(2*order+3))


mean(eqkernel(1,:)')

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
    for i=1:(L+1)
       %%%Have to check why to big...
        if(i==0)
             ha(:,i)= Cdp*sigma'./( mean(Tobs)*( dens.*sum( ( regsx*1 ).^2  )));
        else
             ha(:,i)= Cdp*sigma'*(Vduala(:,i).^2)./( mean(Tobs)*( dens.*sum( ( regsx*Vduala(:,i) ).^2  )));  
        end
           
    end
    H= ha'.^(1/(2*order+2+g) ) ;
end

Xsmoa   =[] ;   
Xsmo2a  =[] ;

for i=0:L
    Ynew=[];
    for j=1:N
        if(i==0)
            Ynew= [Ynew, cell2mat( Y(j) )'*1/N];            
        else
            Ynew= [Ynew, cell2mat( Y(j) )'*Vduala(j,i)];
        end
    end

    [W_dm,fit]=multiloc( cell2mat( Xc' ) , N*Ynew' ,x ,H(i+1),order, kernel  );
     Xsmo2a  =[Xsmo2a fit(:,pos2)];
end


Meansmo2b   =Xsmo2a(:,1) ;                %Mean curve
size(Xsmo2a(:,2:(L+1) ))
V2a         =Xsmo2a(:,2:(L+1) )*sqrt(Da(1:L,1:L))^-1;	%Eigenvectors 2nd derivative smoothed estimate
loadsa      = sqrt(Da)*Vduala';             %Loadings


%%%Final estimator
hX2r        =repmat(Meansmo2b,[1 N])+ V2a(:,1:L)*loadsa(1:L,:);

%hX2r        =repmat(Meansmo2b,[1 N])+ (Xsmo2a-repmat(Meansmo2b,[1 N]) ) *Vduala(:,1:L)*Vduala(:,1:L)';
%hX2r        =repmat(Meansmo2b,[1 N])+ Xa *Vduala(:,1:L)*Vduala(:,1:L)';

end