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




function [Ws] = equivkernel(X,p,kernel_c)

%%%Derive needed Parametes
%p=p-1;
g= size(X,2) ;
T = size(X,1);
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

     Xn = X; 
 
    Z = ones( T, pp );   %%%Construct polynomial
    for m=1:pp
         Z(:,m)= prod( Xn.^repmat(k(m,:),T,1) ,2);
         %Z(:,m*g+1:(m+1)*g) = Xn.^repmat(k(m+1,:),T,1);  
    end
    
    %%Construct Bandwidth matrix
    %H=diag(h);
    
    %%Construct g-dimensional product kernel    
    W= kernel(Xn, kernel_c);  
    W=diag(W) ;
    
A=Z' * W * Z;
%Ws= (-min(X)+max(X))^2* pinv(A) * Z' * W ;
Ws=T/prod(range(X)) * pinv(A) * Z' * W ;
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
