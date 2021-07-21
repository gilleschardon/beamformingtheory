%% Bias, variance, MSE 

% Fig 1 position
% Fig 2 amplitude at x0
% Fig 3 amplitude at a fixed distance

%% Array

%number of microphones
N = 81;

[ya1, za1] = meshgrid(linspace(-0.5, 0.5, sqrt(N)), linspace(-0.5, 0.5, sqrt(N)));

Array = [0*ones(N, 1), ya1(:), za1(:)];

% reference points
X0 = [0,0,0];
X02 = [0.,-0.375,0];

% source
XS = [1.5,0.5,0];

% SNR and frequency
SNR = 10;
KK = 6:2:30;
% Number of realizations
Nt = 1000;

% initialization grid
xg = linspace(0.5, 2.5, 40);
yg = linspace(0, 1, 20);
zg = 0;
[XG YG ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];

% bounds for the optimization
LB = [0.5, 0, -0.5];
UB = [2.5, 1, 0.5];

%% Simulations
[E, B, V,Ep, Bp, Vp,Eps, Bps, Vps] = MCbf_cond(Array, X0, XS, Nt, SNR, KK, Xg, LB, UB);
[E2, B2, V2,Ep2, Bp2, Vp2,Eps2, Bps2, Vps2] = MCbf_cond(Array, X02, XS, Nt, SNR, KK, Xg, LB, UB);

save FIGplanar
%% Position
load FIGplanar


styleI = '-^';
styleII = '-.+';
styleIII = '-.x';
styleIV = '-v';
styleIIb = '--+';
styleIIIb = '--x';

fsize = [100, 100, 500, 800];

MS = 15;
figure

set(gcf, 'Position',  fsize);

He = KK * sqrt(2)/2/pi;

subplot(3, 1, 1)

loglog(He, squeeze(E(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
loglog(He, squeeze(E(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(E(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(E(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(E2(:, :, 2)), styleIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(E2(:, :, 3)), styleIIIb, 'linewidth', 2, 'Markersize', MS)

%ylim([0.00025 5])


ylabel('MSE (m^2)')
xlabel('He')

subplot(3, 1, 2)

loglog(He, squeeze(V(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
loglog(He, squeeze(V(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(V(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(V(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(V2(:, :, 2)), styleIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(V2(:, :, 3)), styleIIIb, 'linewidth', 2, 'Markersize', MS)
ylim([0.001 0.07])
legend('I', 'II', 'III', 'IV', "II'", "III'")

ylabel('Variance (m^2)')
xlabel('He')

subplot(3, 1, 3)

loglog(He, squeeze(B(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
loglog(He, squeeze(B(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(B(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(B(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(B2(:, :, 2)), styleIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(B2(:, :, 3)), styleIIIb, 'linewidth', 2, 'Markersize', MS)
%ylim([1e-8  10])

ylabel('Bias squared (m^2)')
xlabel('He')

%% Amplitude at the reference point
figure
set(gcf, 'Position',  fsize);

He = KK * sqrt(2)/2/pi;

subplot(3, 1, 1)

loglog(He, squeeze(Ep(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
loglog(He, squeeze(Ep(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Ep(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Ep(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+2);

loglog(He, squeeze(Ep(:, :, 5)), ':o', 'linewidth', 2, 'Markersize', MS)


%ylim([0.01, 1.450])

ylabel('MSE (Pa^2)')
xlabel('He')

subplot(3, 1, 2)

loglog(He, squeeze(Vp(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
loglog(He, squeeze(Vp(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Vp(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Vp(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+2);


loglog(He, squeeze(Vp(:, :, 5)), ':o', 'linewidth', 2, 'Markersize', MS)

%ylim([0.02, 0.050])
xlabel('He')

ylabel('Variance (Pa^2)')
legend('I', 'II', 'III', 'IV', "MLE")


subplot(3, 1, 3)

semilogx(He, squeeze(Bp(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
plot(He, squeeze(Bp(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
plot(He, squeeze(Bp(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
plot(He, squeeze(Bp(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+2);

plot(He, squeeze(Bp(:, :, 5)), ':o', 'linewidth', 2, 'Markersize', MS)


ylim([-0.024 0.05])
ylabel('Bias (Pa)')
xlabel('He')

%% Amplitude at 1m
figure
set(gcf, 'Position',  fsize);

He = KK * sqrt(2)/2/pi;

subplot(3, 1, 1)

hold on

%loglog(He, squeeze(Eps(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
loglog(He, squeeze(Eps(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)

loglog(He, squeeze(Eps(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
%loglog(He, squeeze(Eps(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(He, squeeze(Eps2(:, :, 2)), styleIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Eps2(:, :, 3)), styleIIIb, 'linewidth', 2, 'Markersize', MS)

loglog(He, squeeze(Eps(:, :, 5)), ':o', 'linewidth', 2, 'Markersize', MS)

%ylim([0.15, 1])
xlabel('He')

ylabel('MSE (Pa^2)')

subplot(3, 1, 2)

%loglog(He, squeeze(Vps(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
%hold on
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

loglog(He, squeeze(Vps(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Vps(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
%loglog(He, squeeze(Vps(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(He, squeeze(Vps2(:, :, 2)), styleIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Vps2(:, :, 3)), styleIIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Vps(:, :, 5)), ':o', 'linewidth', 2, 'Markersize', MS)

%ylim([0.08, 0.55])
xlabel('He')

ylabel('Variance (Pa^2)')
legend('II', 'III', "II'", "III'", "MLE")

subplot(3, 1, 3)

%plot(He, squeeze(Bps(:, :, 1)), styleI, 'linewidth', 2, 'Markersize', MS)
hold on
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

set(gca, 'XScale', 'log')
semilogx(He, squeeze(Bps(:, :, 2)), styleII, 'linewidth', 2, 'Markersize', MS)
hold on
plot(He, squeeze(Bps(:, :, 3)), styleIII, 'linewidth', 2, 'Markersize', MS)
%plot(He, squeeze(Bps(:, :, 4)), styleIV, 'linewidth', 2, 'Markersize', MS)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(He, squeeze(Bps2(:, :, 2)), styleIIb, 'linewidth', 2, 'Markersize', MS)
loglog(He, squeeze(Bps2(:, :, 3)), styleIIIb, 'linewidth', 2, 'Markersize', MS)
plot(He, squeeze(Bps(:, :, 5)), ':o', 'linewidth', 2, 'Markersize', MS)


%ylim([-1, 1.1])

ylabel('Bias (Pa)')
xlabel('He')
