function [B, dref] = objB3cond(X, Xm, X0, sig, k)

 dref = norm(X-X0);

 ax = dref * dictionary(Xm, X, k);
 bb = ax./norm(ax)^2;
 
 B = abs(bb'*sig)^2;

 return