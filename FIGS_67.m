
%% author Gilles Chardon, 2021

%% Bias of formulation II


% Array
M1 = 25;

[ya1, za1] = meshgrid(linspace(-0.5, 0.5, sqrt(M1)), linspace(-0.5, 0.5, sqrt(M1)));

% z coordinate of the unbiased source
XS0 = 4;

% z coordinates of the sources
XX = [2 3.5 4 5 6 7];

% z coordinates of the microphones
xx0 = (XS0 - sqrt(XS0^2 - 4 * (ya1(:).^2 + za1(:).^2)))/2;

Array = [xx0 , ya1(:), za1(:)];
    
%% Geometry, array and sources
figure

ya1 = linspace(-0.5, 0.5, sqrt(M1))';
xxa0 = (XS0 - sqrt(XS0^2 - 4 * (ya1(:).^2)))/2;

scatter(xxa0, ya1, 10, 'ok', 'filled')
hold on
scatter(XX, XX*0, 'xk')
scatter(XS0, 0, 200, 'xk', 'linewidth', 3)
axis image

xlim([-0.2, 7.2])
ylim([-0.7, 0.7])
xlabel('z (m)')
ylabel('y (m)')

%%
    
X00 = [0,0,0];
 
SNR = 0;
KK = 15:1:50;

He = KK * sqrt(2)/2/pi;



Z1= zeros(length(KK), 2);
Z2 = zeros(length(KK), 2);
Z3 = zeros(length(KK), 2);
Z4 = zeros(length(KK), 2);

bias1 = zeros(length(KK), length(XX));
bias2 = zeros(length(KK), length(XX));
bias3 = zeros(length(KK), length(XX));
bias4 = zeros(length(KK), length(XX));
biasest = zeros(length(KK), length(XX));


for narray = 1:length(XX)
    

for nk = 1:length(KK)

    XS = [XX(narray),0,0];

    
    X0 =X00;
    k = KK(nk);

a = dictionary(Array, XS, k);

a0 = dictionary(X0, XS, k);

sigs0 = a;


sigs = sigs0;

funB1  = @(x) - objB1cond([x 0], Array, X0, sigs, k);
funB2  = @(x) - objB2cond([x 0], Array, X0, sigs, k);
funB3  = @(x) - objB3cond([x 0], Array, X0, sigs, k);
funB4  = @(x) - objB4cond([x 0], Array, X0, sigs, k);

[Z1(nk, :), P1] = fminunc(funB1, XS(1:2));
[Z2(nk, :), P2] = fminunc(funB2, XS(1:2));
[Z3(nk, :), P3] = fminunc(funB3, XS(1:2));
[Z4(nk, :), P4] = fminunc(funB4, XS(1:2));


% approximated bias
Deltai = Array(:, 1) - XX(narray);
Delta0 =  - XX(narray);
Ri = sqrt( (Array(:, 1) - XX(narray)).^2 + Array(:, 2).^2 + Array(:, 3).^2);
R0 = XX(narray);

dd = - 2 * sum (Deltai./Ri.^2 - Delta0./R0.^2) / M1;
dd2 =  (sum(Deltai./Ri).^2/M1 - sum((Deltai./Ri).^2)) * k^2 / M1 * 2; 

biasest(nk, narray) = - dd/dd2;

end

bias1(:, narray) = Z1(:, 1) - XX(narray);
bias2(:, narray) = Z2(:, 1) - XX(narray);
bias3(:, narray) = Z3(:, 1) - XX(narray);
bias4(:, narray) = Z4(:, 1) - XX(narray);


end
%%
figure
hold on
plot(He, fliplr(bias2), 'linewidth', 2)
set(gca,'ColorOrderIndex',1)

plot(He, fliplr(biasest))

xlabel('He')
ylabel('bias (m)')

xlim([min(He) max(He)])
ylim([-0.3 0.5])


legend(compose("z_s = %.1f", fliplr(XX)))
