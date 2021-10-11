
%% author Gilles Chardon, 2021

clear all
M = 3;


[xa, ya, za] = meshgrid(linspace(-0.5, 0.5, M), linspace(-0.5, 0.5, M), linspace(-0.5, 0.5, M));

Array = [xa(:), ya(:), za(:)];



% removing the center
Array(14, :) = [];
X0 = [0,0,0];

XS = [0.1 0.05 0.45];

SNR = 10;
KK = 2:0.5:10;

Nt = 1000;


xg = linspace(-0.5, 0.5, 20);
yg = linspace(-0.5, 0.5, 20);
zg = linspace(-0.5, 5, 20);
[XG YG ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];

%%
L = length(KK);
zz = linspace(0, 1.5, L);
LB = [-0.5 -0.5 -0.5];
UB = [0.5 0.5 0.5];
for u = 1:L
XS = [0.1 0.3 0.4];

[E, B, V,Ep, Bp, Vp,Eps, Bps, Vps] = MCbf_cond(Array, X0, XS, Nt, SNR, KK(u), Xg, LB, UB);

E1(u) = E(1);
E2(u) = E(2);
E3(u) = E(3);
E4(u) = E(4);
end

save FIGcubic
%%
close all
load FIGcubic

styleI = '-^';
styleII = '-.+';
styleIII = '-.x';
styleIV = '-v';
styleIIb = '--+';
styleIIIb = '--x';


semilogy(KK, squeeze(E1), styleI', 'linewidth', 2, 'Markersize', 10)
hold on
loglog(KK, squeeze(E2), styleII, 'linewidth', 2, 'Markersize', 10)
loglog(KK, squeeze(E3), styleIII, 'linewidth', 2, 'Markersize', 10)
loglog(KK, squeeze(E4), styleIV, 'linewidth', 2, 'Markersize', 10)

legend('I', 'II', 'III', 'IV')
ylabel('MSE (m^2)')
xlabel('wavenumber k (m^{-1})')

figure

plot(E1./E4)
