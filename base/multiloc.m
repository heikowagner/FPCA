% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: g-Dimensional local polynomial estimator 
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
%X - g dimensional inputcoordinate , 
%Y  - function value
%x - g dimensional output coordinate
%h - gxg bandwidth matrix dimensional vector of  bandwiths 
%p  - degree of polynomial
 
% ------------------------------------------------------------------------------ 
% Output:      
%fit  - Smoothed curves and derivatives
%W    - Weight used for diagonal correction

% ------------------------------------------------------------------------------ 
% Keywords:    local polynomial surface estimator, derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2017/01/24
% ------------------------------------------------------------------------------ 




function [W,fit] = multiloc( X,  Y, x, H ,p,kernel_c)

%%%Derive needed Parametes
%p=p-1;
g= size(X,2) ;
T = size(X,1);
T_out = size(x,1);
dummy=permn(0:p,g);             %%%Get all possible combinations
k=dummy( sum(dummy,2)<p+1,: );
pp=size(k,1);

function [A]=kernel(Xn, k_type)
g=size(Xn,2);
A=1;
    if(strcmp(k_type,'Gauss')==1)   %Use Gaussian Kernel
        for (m=1:g)
           A= A.*normpdf(Xn(:,m))  ;
        end
    else                            %Use Epanechnikov Kernel
        for (m=1:g)
           A= A.*epan(Xn(:,m))  ;
        end
    end   
end

function [out]=getpoint(x_j)
    Xn=zeros( size(X) );
    Xn = X -  repmat( x_j,T,1 ); 
 
    Z = ones( T, pp );   %%%Construct polynomial
    for m=1:pp
         Z(:,m)= prod( Xn.^repmat(k(m,:),T,1) ,2);
         %Z(:,m*g+1:(m+1)*g) = Xn.^repmat(k(m+1,:),T,1);  
    end
    
    %%Construct Bandwidth matrix
    %H=diag(h);
    
    %%Construct g-dimensional product kernel    
    W= det(H)^(-1)*kernel(Xn*H^(-1), kernel_c);  
    W=diag(W) ;
    
A=Z' * W * Z;
%Ws= pinv(A, 0.0001) * Z' * W;
A;
Ws= inv(A) * Z' * W;
b = Ws * Y;

fit_j = ones( 1, pp );
    for m=1:pp
        fit_j(:,m ) = b(m,1);
    end

W_j = ones( 1, pp );
     for m=1:pp
        W_j(:,m ) =  sum(Ws(m,:).^2);
    end
 
out=[W_j,fit_j]';
end 
[out]=arrayfun(@(s) getpoint( x(s,:) ), 1:T_out , 'UniformOutput', false);
dummy=cell2mat(out)';
scaler = max(1,factorial( sum(k,2) ));
W=dummy(:,1:pp)*diag(scaler.^2);
fit=dummy(:,(pp+1):end)*diag(scaler);
end


% 
% %%tryouts
% p=3
% g=2
% 
% %%get all possible combinations
% 
% v=repmat(0:p,1,g)
% perms(1:p)
% 
% v'v
% 
% combnk(0:3,2)
% 
% perms(0:2)
% 
% X=repmat([2,3],20,1)
% d=[1,2]
% 
% X.^repmat(d,20,1)   %%%Das klappt schonmal... schön^^ 
% 
% 
