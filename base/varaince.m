% ------------------------------------------------------------------------------ 
% Article:     Functional Principal Component Analysis for Derivatives of 
%              High-Dimensional Spatial Curves
% ------------------------------------------------------------------------------ 
% Description: 2 Dimensional Variance estimator, Based on Hall, Marron 1990
%
% ------------------------------------------------------------------------------ 
% Usage:       - 
% ------------------------------------------------------------------------------ 
% Inputs Simulation:      
%X - g-dimesnional spacial coordinates, 
%Y  - function value
 
% ------------------------------------------------------------------------------ 
% Output:      
%sigma

% ------------------------------------------------------------------------------ 
% Keywords:    local polynomial surface estimator, derivatives
% ------------------------------------------------------------------------------ 
% See also:    -  
% ------------------------------------------------------------------------------ 
% Author:      Heiko Wagner, 2016/21/21 
% ------------------------------------------------------------------------------ 


function [sigma ] = variance( X, Y)

 %i=1
 %X1=cell2mat(  c5unil(i) )
 %X2= cell2mat(  c2unil(i) ) 
 %Y=cell2mat(  c3unil(i) )
 
g= size(X,2);

kernel='';

N1 = size(X,1);
T = N1;

h0=T^(-2/(4+g));
 
Wges=zeros(N1,T);

for (j=1:T)
    Xmx=zeros( size(X) );
    for (m=1:g) %%%construct g dimensional grid around t_i
        Xmx(:,m) = X(:,m) - X(j,m);
    end
    if(strcmp(kernel,'Gauss')==1)   
        %W= normpdf(Xmx1/h0) .*normpdf(Xmx2/h0)  ;
        W=1;
        for (m=1:g) %%%construct g dimensional grid around t_i
           W= W.*normpdf(Xmx(:,m)/h0)  ;
        end
        %'Using Gaussian Kernel'
    else
        W=1;
        for (m=1:g) %%%construct g dimensional grid around t_i
           W= W.*epan(Xmx(:,m)/h0)  ;
        end        
    end
Wges(:,j)=W/sum(W);

end
v= T - 2 * sum(diag(Wges)) + sum(sum( Wges.^2)) ;
sigma= 1/v* sum( (Y-Wges'*Y).^2 );
end

%
% scatter3(X1,X2,Ws(1,:))
% 
% scatter3(my,mx,fit0)
% hold on
% scatter3(X1,X2,Y)
% 
% scatter(X2,epan(Xmx2/h2))
% hold on
% scatter(X1,epan(Xmx1/h1))
% 
% 
% x1(j)
% x2(j)

%varaince(cell2mat(  c5unil(i) ),  cell2mat(  c2unil(i) ) ,  cell2mat(  c3unil(i) ),my,mx,kernel='epan') 

%varaince(cell2mat(  c5unil(i) ),  cell2mat(  c2unil(i) ) ,  cell2mat(  c3unil(i) ),my,mx) 
