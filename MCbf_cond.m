function [E, B, V,Ep, Bp, Vp,Eps, Bps, Vps] = MCbf_cond(Array, X0, XS, Nt, SNR, KK, Xg, LB, UB)

  %% author Gilles Chardon, 2021
  
Nmethods = 4;

Nk = length(KK);
Nsnr = length(SNR);

% MSE, Variance, Bias
% Position, Amplitude at the reference point, at 1m
E = zeros(Nk, Nsnr, Nmethods);
B = zeros(Nk, Nsnr, Nmethods);
V = zeros(Nk, Nsnr, Nmethods);
Ep = zeros(Nk, Nsnr, Nmethods+1);
Bp = zeros(Nk, Nsnr, Nmethods+1);
Vp = zeros(Nk, Nsnr, Nmethods+1);
Eps = zeros(Nk, Nsnr, Nmethods+1);
Bps = zeros(Nk, Nsnr, Nmethods+1);
Vps = zeros(Nk, Nsnr, Nmethods+1);

NN = Nk*Nsnr*Nt;

m = 0;

Astar = 1;
Aref = Astar / norm(XS-X0);

for nk = 1:length(KK)

    
    k = KK(nk);

    a = dictionary(Array, XS, k)*Astar;

    a0 = dictionary(X0, XS, k);

    
    Dg = dictionary(Array, Xg, k)*norm(XS-X0);
    
    % steering vectors of the grid
    DI = Dg./abs(Dg);
    DII = Dg./abs(Dg.^2);
    DIII = Dg./sum((abs(Dg).^2), 1);
    DIV = Dg./sqrt(sum((abs(Dg).^2), 1));


Z1 = zeros(Nt, 3);
Z2 = zeros(Nt, 3);
Z3 = zeros(Nt, 3);
Z4 = zeros(Nt, 3);


    P0 = abs(a0)^2;

    for nsnr = 1:Nsnr
        for nt = 1:Nt
            
            waitbar(m/NN)
            m = m + 1;


            % data generation
            sigs0 = a;

            SNR0 = SNR(nsnr);
    
            noise = randn(size(sigs0)) + 1i*randn(size(sigs0));
            noise = noise / norm(noise, 'fro') * norm(sigs0, 'fro') * 10^(-SNR0/20);

            sigs = sigs0 + noise;
        
            funB1  = @(x) - objB1cond(x, Array, X0, sigs, k);        
            funB2  = @(x) - objB2cond(x, Array, X0, sigs, k);
            funB3  = @(x) - objB3cond(x, Array, X0, sigs, k);
            funB4  = @(x) - objB4cond(x, Array, X0, sigs, k);
            
            % estimations
            
            [~, idx] = max( abs(DI'*sigs0));
            [Z1(nt, :)] = fmincon(funB1, Xg(idx, :), [], [], [], [], LB, UB);
            [o, d] = objB1cond(Z1(nt, :), Array, X0, sigs, k);
            P1(nt) = sqrt(o);
            P1s(nt) = P1(nt) * d;

            [~, idx] = max( abs(DII'*sigs0));
            [Z2(nt, :)] = fmincon(funB2, Xg(idx, :), [], [], [], [], LB, UB);
            [o, d] = objB2cond(Z2(nt, :), Array, X0, sigs, k);
            P2(nt) = sqrt(o);
            P2s(nt) = P1(nt) * d;

            [~, idx] = max( abs(DIII'*sigs0));
            [Z3(nt, :)] = fmincon(funB3, Xg(idx, :), [], [], [], [], LB, UB);
            [o, d] = objB3cond(Z3(nt, :), Array, X0, sigs, k);
            P3(nt) = sqrt(o);
            P3s(nt) = P3(nt) * d;

            [~, idx] = max( abs(DIV'*sigs0));
            [Z4(nt, :)] = fmincon(funB4, Xg(idx, :), [], [], [], [], LB, UB);
            [o, d] = objB4cond(Z4(nt, :), Array, X0, sigs, k);
            P4(nt) = sqrt(o);
            P4s(nt) = P4(nt) * d;

            [o, d] = objB3cond(Z4(nt, :), Array, X0, sigs, k);
            P5(nt) = sqrt(o);
            P5s(nt) = P5(nt) * d;        

        end

    E(nk, nsnr,  1) = norm(Z1 - XS(1:3), 'fro')^2 / Nt;
    E(nk, nsnr,  2) = norm(Z2 - XS(1:3), 'fro')^2 / Nt;
    E(nk, nsnr,  3) = norm(Z3 - XS(1:3), 'fro')^2 / Nt;
    E(nk, nsnr,  4) = norm(Z4 - XS(1:3), 'fro')^2 / Nt;

    B(nk, nsnr, 1) = norm(mean(Z1, 1) - XS(1:3))^2;
    B(nk, nsnr,  2) = norm(mean(Z2, 1) - XS(1:3))^2;
    B(nk, nsnr,  3) = norm(mean(Z3, 1) - XS(1:3))^2;
    B(nk, nsnr,  4) = norm(mean(Z4, 1) - XS(1:3))^2;

    V(nk, nsnr, 1) = sum(var(Z1, 1, 1));
    V(nk, nsnr,  2) = sum(var(Z2, 1, 1));
    V(nk, nsnr,  3) = sum(var(Z3, 1, 1));
    V(nk, nsnr, 4) = sum(var(Z4, 1, 1));
   
    Ep(nk, nsnr,  1) = norm(P1 - Aref)^2 / Nt;
    Ep(nk, nsnr,  2) = norm(P2 - Aref)^2 / Nt;
    Ep(nk, nsnr,  3) = norm(P3 - Aref)^2 / Nt;
    Ep(nk, nsnr,  4) = norm(P4 - Aref)^2 / Nt;
    Ep(nk, nsnr,  5) = norm(P5 - Aref)^2 / Nt;

    Bp(nk, nsnr, 1) = (mean(P1) - Aref);
    Bp(nk, nsnr,  2) = (mean(P2) - Aref);
    Bp(nk, nsnr,  3) = (mean(P3) - Aref);
    Bp(nk, nsnr,  4) = (mean(P4) - Aref);
    Bp(nk, nsnr,  5) = (mean(P5) - Aref);

    Vp(nk, nsnr, 1) = (var(P1));
    Vp(nk, nsnr,  2) =(var(P2));
    Vp(nk, nsnr,  3) =(var(P3));
    Vp(nk, nsnr, 4) = (var(P4));
    Vp(nk, nsnr, 5) = (var(P5));
        
    Eps(nk, nsnr,  1) = norm(P1s - Astar)^2 / Nt;
    Eps(nk, nsnr,  2) = norm(P2s - Astar)^2 / Nt;
    Eps(nk, nsnr,  3) = norm(P3s - Astar)^2 / Nt;
    Eps(nk, nsnr,  4) = norm(P4s - Astar)^2 / Nt;
    Eps(nk, nsnr,  5) = norm(P5s - Astar)^2 / Nt;

    Bps(nk, nsnr, 1) = (mean(P1s) - Astar);
    Bps(nk, nsnr,  2) = (mean(P2s) - Astar);
    Bps(nk, nsnr,  3) = (mean(P3s) - Astar);
    Bps(nk, nsnr,  4) = (mean(P4s) - Astar);
    Bps(nk, nsnr,  5) = (mean(P5s) - Astar);

    Vps(nk, nsnr, 1) = (var(P1s));
    Vps(nk, nsnr,  2) =(var(P2s));
    Vps(nk, nsnr,  3) =(var(P3s));
    Vps(nk, nsnr, 4) = (var(P4s));
    Vps(nk, nsnr, 5) = (var(P5s));

    end
end
end

