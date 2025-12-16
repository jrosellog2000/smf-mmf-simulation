lambda = 1.55e-6;       % wavelength in meters
a = 4.1e-6;             % core radius in meters
n_core = 1.450;         
n_clad = 1.444;         

NA = sqrt(n_core^2 - n_clad^2);   % numerical aperture
V = 2*pi*a/lambda*NA;             % normalized frequency

N = 512;        % number of grid points
L = 7e-6;      % grid size in meters
x = linspace(-L,L,N);
y = x;
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);

w0 = a/sqrt(2);            % mode field radius
E_smf = exp(-(R.^2)/(w0^2));
E_smf = E_smf/max(E_smf(:));

% Plot intensity
figure;
imagesc(x*1e6,y*1e6,abs(E_smf).^2)
axis equal tight
xlabel('µm'); ylabel('µm');
title('SMF LP_{01} Mode');
colorbar;
