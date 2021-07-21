function MCbf_line(Array, X0, XS, k, yy, N, X02)

a = dictionary(Array, XS, k) * norm(XS-X0);
zz = 0;
xx = XS(1);
[Xg Yg Zg] = meshgrid(xx, yy, zz);
Xg = [Xg(:), Yg(:), Zg(:)];



[B1, B2, B3, B4] = BF1234_cond(a, Array, Xg, X0, k);
[B12, B22, B32, B42] = BF1234_cond(a, Array, Xg, X02, k);

               
hh = hot;
hh = hh(1:220, :);

styleI = '-';
styleII = '-.';
styleIII = '-.';
styleIV = '-';
styleIIb = '--';
styleIIIb = '--';

figure
subplot(1, 3, 1)
hold on

plot(yy, B1,styleI, 'linewidth', 2)
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+2);

plot(yy, B4, styleIV, 'linewidth', 2)%, 'linecolor', 'none')
ylabel('B_U')

[~, idx] = max(B1);
scatter(yy(idx), B1(idx), 'k', 'filled');
[~, idx] = max(B4);
scatter(yy(idx), B4(idx), 'k', 'filled');
grid on
legend('I', 'IV', 'location', 'southeast')

xlabel('y (m)')

subplot(1, 3, 2)
hold on

set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

plot(yy, B2, styleII, 'linewidth', 2)%, 'linecolor', 'none')

plot(yy, B3, styleIII, 'linewidth', 1)%, 'linecolor', 'none')


[~, idx] = max(B2);
scatter(yy(idx), B2(idx), 'k', 'filled');

[~, idx] = max(B3);
scatter(yy(idx), B3(idx), 'k', 'filled');

legend('II', 'III', 'location', 'southeast')

grid on
xlabel('y (m)')
ylabel('B_U')

subplot(1, 3, 3)
hold on

set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);
plot(yy, B22, styleIIb, 'linewidth',2)%, 'linecolor', 'none')

plot(yy, B32, styleIIIb', 'linewidth', 1)%, 'linecolor', 'none')


ylabel('B_U')

[~, idx] = max(B22);
scatter(yy(idx), B22(idx), 'k', 'filled');

[~, idx] = max(B32);
scatter(yy(idx), B32(idx), 'k', 'filled');

xlabel('y (m)')
legend("II'", "III'", 'location', 'southeast')

grid on
