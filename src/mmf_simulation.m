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
[X,Y] = meshgrid(x,y); % obtaining each coordinate
R = sqrt(X.^2 + Y.^2); % distance from the centre to each point
phi = atan2(Y, X); % azimuthal angle for each coordinate


% Using E_lp_modes function to obtain the complex electrical field for each LP mode
E01 = E_lp_modes(0,1,phi,R,a);
E11 = E_lp_modes(1,1,phi,R,a);
E21 = E_lp_modes(2,1,phi,R,a);
E02 = E_lp_modes(0,2,phi,R,a);
E31 = E_lp_modes(3,1,phi,R,a);
E12 = E_lp_modes(1,2,phi,R,a);
E41 = E_lp_modes(4,1,phi,R,a);


% Combine modes
r = rand(1,7);
E_mmf = r(1)*E01 + r(2)*E11 + r(3)*E21 + r(4)*E31 + r(5)*E12 + r(6)*E02 + r(7)*E41;
E_mmf = E_mmf/max(abs(E_mmf(:)));

% Plot intensity
figure;
imagesc(x*1e6,y*1e6,abs(E_mmf).^2) % intensity os is E^2
axis equal tight
xlabel('µm'); ylabel('µm');
title('MMF LP Modes');
colorbar;
