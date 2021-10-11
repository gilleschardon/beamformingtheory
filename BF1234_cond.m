function [B1, B2, B3, B4] = BF1234_cond(sig, Xm, Xs, X0, k)

  %% author Gilles Chardon, 2021
  
D = dictionary(Xm, Xs, k);

dref = sqrt((Xs(:, 1)-X0(1)).^2 + (Xs(:, 2)-X0(2)).^2 + (Xs(:, 3)-X0(3)).^2)';
Dref = D .* dref;

normsDref2 = sum(abs(Dref).^2, 1);

N = size(D, 1);

Bm1 = 1/N * Dref ./ abs(Dref);
Bm2 = 1/N * Dref ./ abs(Dref).^2;
Bm3 = Dref ./ normsDref2;
Bm4 = 1/sqrt(N) * Dref ./ sqrt(normsDref2);


B1 = abs(Bm1'*sig).^2;
B2 = abs(Bm2'*sig).^2;
B3 = abs(Bm3'*sig).^2;
B4 = abs(Bm4'*sig).^2;



end
