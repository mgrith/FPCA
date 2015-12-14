function [k]=epan(t);
k= 3/4*(1-t.^2);
k((abs(t)-1)>0)=0;
end