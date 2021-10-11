function [B, dref] = objB3cond(X, Xm, X0, sig, k)

  %% author Gilles Chardon, 2021
  
 dref = norm(X-X0);

 ax = dref * dictionary(Xm, X, k);
 bb = ax./norm(ax)^2;
 
 B = abs(bb'*sig)^2;

 return
