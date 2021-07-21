
%% Bias in 2D and 1D

% Figures 1 2 and 3

%% Array and grid

% number of microphones
N = 81;

% coordinates
[ya1, za1] = meshgrid(linspace(-0.5, 0.5, sqrt(N)), linspace(-0.5, 0.5, sqrt(N)));

Array = [0*ones(N, 1), ya1(:), za1(:)];

% source
XS = [1.5,0.5,0];

% SNR (noiseless measurements)
SNR = 10000;

k2D = 15;

% size of the 2D grid
N = 101;

yy = linspace(0, 1, N);
xx = linspace(0.5, 2.5, 2*N);

% reference point 1
X0 = [0,0,0];
% reference point 2
X02 = [0.,-0.375,0];

%% Fig 1. geometry
close all
figure

scatter(Array(:, 1), Array(:, 2), 'ok', 'filled')
hold on
rectangle('Position', [0.5 0 2 1], 'FaceColor',[0.8 .8 .8], 'EdgeColor', [1 1 1])
line([1.5 1.5], [0.3, 0.6], 'Color', 'black')
scatter(X0(1), X0(2), 100, 'sk')
scatter(X02(1), X02(2), 100, '^k')
scatter(XS(1), XS(2), 100, '*k')


axis equal
ylim([-0.6, 1.1])
xlim([-0.1, 2.6])
xlabel('z (m)')
ylabel('y (m)')
legend('microphones', 'x_0', "x_0'", 'Source', 'location', 'southeast')


%% Fig 2. Contour plot, two references
MCbf_contour(Array, X0, XS, k2D, xx, yy, N, 2*N, X02);

%% Fig 3/ Criteria on a line

k1D = 5;

yy2 =  linspace(0.3, 0.6, N);

MCbf_line(Array, X0, XS, k1D, yy2, N, X02);
