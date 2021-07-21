function [B, dref] = objB1cond(X, Xm, X0, sig, k)

 dref = norm(X-X0);

 ax = dref * dictionary(Xm, X, k);
 bb = 1/size(Xm, 1) * ax./abs(ax);
 
 B = abs(bb'*sig)^2;

 return