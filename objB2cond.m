function [B, dref] = objB2cond(X, Xm, X0, sig, k)

 dref = norm(X-X0);

 ax = dref * dictionary(Xm, X, k);
 bb = 1/size(Xm, 1) * ax./abs(ax).^2;
 
 B = abs(bb'*sig)^2;

 return