function [mse] = mse_ord( X,  Y, p)

%mse=mean(( sign(X(p))*X - sign(Y(p))*Y).^2 ) 

mse=min( mean(( X - Y).^2 ) ,mean(( X + Y).^2 ))

end
