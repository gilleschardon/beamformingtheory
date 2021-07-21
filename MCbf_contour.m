function [B1, B2, B3, B4] = MCbf_contour(Array, X0, XS, k, xx, yy, N, M, X02)

% source
a = dictionary(Array, XS, k) * norm(XS-X0);
zz = 0;

%vgrid
[Xg Yg Zg] = meshgrid(xx, yy, zz);

Xg = [Xg(:), Yg(:), Zg(:)];


% values of the criteria on the grid
[B1, B2, B3, B4] = BF1234_cond(a, Array, Xg, X0, k);
[B12, B22, B32, B42] = BF1234_cond(a, Array, Xg, X02, k);


%% figures

nlev = 13;
figure
subplot(3,2,1)

nlev = linspace(0, 1.2, nlev);
hh = hot;
hh = hh(1:220, :);


contourf(xx, yy, reshape(B1, N,M), nlev)%, 'linecolor', 'none')
hold on
[~, idx] = max(B1);
scatter(Xg(idx, 1), Xg(idx, 2), 'k', 'filled')
scatter(XS(1), XS(2), 100, 'ok', 'linewidth', 2);
colormap(hh)
colorbar
axis image
title('I')
xlabel('z (m)')
ylabel('y (m)')



subplot(3,2,3)
contourf(xx, yy, reshape(B2, N,M),nlev)%, 'linecolor', 'none')
hold on
[~, idx] = max(B2);
scatter(Xg(idx, 1), Xg(idx, 2), 'k', 'filled')
scatter(XS(1), XS(2), 100, 'ok','linewidth', 2);
colormap(hh)
colorbar
axis image

title('II')
xlabel('z (m)')
ylabel('y (m)')

subplot(3,2,4)
contourf(xx, yy, reshape(B3, N,M),nlev)%, 'linecolor', 'none')
hold on
[~, idx] = max(B3);
scatter(Xg(idx, 1), Xg(idx, 2), 'k', 'filled')
scatter(XS(1), XS(2), 100, 'ok', 'linewidth', 2);
colormap(hh)
colorbar
axis image

title("III")
xlabel('z (m)')
ylabel('y (m)')

subplot(3,2,2)
contourf(xx, yy, reshape(B4, N,M),nlev)%, 'linecolor', 'none')
hold on
[~, idx] = max(B4);
scatter(Xg(idx, 1), Xg(idx, 2), 'k', 'filled')
scatter(XS(1), XS(2), 100, 'ok', 'linewidth', 2);
colormap(hh)
colorbar
axis image
title('IV')
xlabel('z (m)')
ylabel('y (m)')

subplot(3,2,5)
contourf(xx, yy, reshape(B22, N,M),nlev)%, 'linecolor', 'none')
hold on
[~, idx] = max(B22);
scatter(Xg(idx, 1), Xg(idx, 2), 'k', 'filled')
scatter(XS(1), XS(2), 100, 'ok', 'linewidth', 2);
colormap(hh)
colorbar
axis image

title("II'")
xlabel('z (m)')
ylabel('y (m)')

subplot(3,2,6)
contourf(xx, yy, reshape(B32, N,M),nlev)%, 'linecolor', 'none')
hold on
[~, idx] = max(B32);
scatter(Xg(idx, 1), Xg(idx, 2), 'k', 'filled')
scatter(XS(1), XS(2), 100, 'ok', 'linewidth', 2);
colormap(hh)
colorbar
axis image


title("III'")
xlabel('z (m)')
ylabel('y (m)')




end
