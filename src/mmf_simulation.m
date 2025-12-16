lambda = 1.55e-6;       % wavelength in meters
a = 50e-6;             % core radius in meters
n_core = 1.450;         
n_clad = 1.444;         

NA = sqrt(n_core^2 - n_clad^2);   % numerical aperture
V = 2*pi*a/lambda*NA;            % normalized frequency
disp(['V-number = ', num2str(V)])

N = 512;        % number of grid points
L = 100e-6;      % grid size in meters
x = linspace(-L,L,N);
y = x;
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);

% Example LP modes
u01 = 2.405;  % LP01 first zero of J01
u11 = 3.832;  % LP11 first zero of J11
u21 = 5.136;  % LP21 first zero of J21
u02 = 5.520;  % LP02 first zero of J02

E01 = besselj(0,u01*R/a).*(R<=a);
E11 = besselj(1,u11*R/a).*cos(atan2(Y,X)).*(R<=a);
E21 = besselj(2,u21*R/a).*cos(2*atan2(Y,X)).*(R<=a);
E02 = besselj(0,u02*R/a).*(R<=a);

% Combine modes
E_mmf = E01 + 0.1*E11 + 0.5*E21 + 0.3*E02;
E_mmf = E_mmf/max(abs(E_mmf(:)));

% Plot intensity
figure;
imagesc(x*1e6,y*1e6,abs(E_mmf).^2)
axis equal tight
xlabel('µm'); ylabel('µm');
title('MMF LP Modes');
colorbar;
