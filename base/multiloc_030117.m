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
%h - g dimensional vector of  bandwiths 
%p  - degree of polynomial
 
% ------------------------------------------------------------------------------ 
% Output:      
%fit0  - Smoothed curves
%fit1  - Smoothed 1. derivative
%fit2  - Smoothed 2. derivative
%fit3  - Smoothed 3. derivative
%h1    - g dimensional vector of  bandwiths 

% ------------------------------------------------------------------------------ 
% Keywords:    local polynomial surface estimator, derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2016/12/22 
% ------------------------------------------------------------------------------ 


function [fit0, fit1, fit2, fit3, W0, W1, W2, W3] = multiloc( X,  Y, x, h,p,kernel)

g=size(X,2);
N1 = size(X,1);
Na = size(x,1);

function [fit0, fit1, fit2, fit3, W0, W1, W2, W3]=getpoint(j)
  
    Xn=zeros( size(X) );
    for (m=1:g) %%%construct g dimensional grid around point j
        Xn(:,m) = X(:,m) - x(j,m);
    end
    
    Z = ones( N1, 2*p+1 );
    for m=1:p
        Z(:,(2*m):(2*m+1) ) = Xn.^m;
    end
    
    %%Construct g-dimensional product kernel
 
    if(strcmp(kernel,'Gauss')==1)   
        %W= normpdf(Xmx1/h0) .*normpdf(Xmx2/h0)  ;
        W=1;
        for (m=1:g) %%%construct g dimensional grid around t_i
           W= W.*normpdf(Xn(:,m)/h(m))/h(m)  ;
        end
        W=diag(W);
        %'Using Gaussian Kernel'
    else
        W=1;
        for (m=1:g) %%%construct g dimensional grid around t_i
           W= W.*epan(Xn(:,m)/h(m))/h(m)  ;
        end 
        W=diag(W);
    end
    

A=Z' * W * Z;
Ws= pinv(A) * Z' * W;
b = Ws * Y;

  fit1 =0;
  W1 = 0;
  fit2 = 0;
  W2 = 0;
  fit3 = 0;
  W3 = 0;
  
%Ws2=diag(Ws*Ws');
 % deriv=0;
  fit0 =   b(1,1);
  W0 =   sum(Ws(1,:).^2) ;
  if(p>1)
  %deriv=1;
  fit1 =factorial(1)*b(3,1); 
  W1 = factorial(1)^2*sum(Ws(3,:).^2);
  end
  if(p>2)
 % deriv=2;
  fit2 =factorial(2) * b(5,1);  
  W2 = factorial(2)^2*sum(Ws(5,:).^2);
  %deriv=3;
  fit3 =factorial(3)*b(7,1); 
  W3 = factorial(2)^2*sum(Ws(7,:).^2);
  end
end
  
[fit0, fit1, fit2, fit3, W0, W1, W2, W3]=arrayfun(@getpoint,1:Na); 

end