function [k]=epan(t);
%t=(-100:100)/50
k= 3/4*(1-t.^2);
   k((abs(t)-1)>0)=0;
end