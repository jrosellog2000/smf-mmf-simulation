
lambda = 1.55e-6;       % wavelength in meters
a = 20e-6;              % core radius in meters
n_core = 1.450;         
n_clad = 1.444;         

NA = sqrt(n_core^2 - n_clad^2);   % numerical aperture
V = 2*pi*a/lambda*NA;             % normalized frequency
disp(['V-number = ', num2str(V)])

N = 512;        % number of grid points
L = 50e-6;      % grid size in meters
x = linspace(-L,L,N);
y = x;
[X,Y] = meshgrid(x,y); 
R = sqrt(X.^2 + Y.^2); 
phi = atan2(Y, X); 

% Matrix of Bessel zeros (rows: l, cols: m)
U = [
    2.4048, 5.5201, 8.6537, 11.7915, 14.9309;
    3.8317, 7.0156, 10.1735, 13.3237, 16.4706;
    5.1356, 8.4172, 11.6198, 14.7960, 17.9598;
    6.3802, 9.7610, 13.0152, 16.2235, 19.4094;
    7.5883, 11.0647, 14.3725, 17.6160, 20.8269;
    8.7715, 12.3386, 15.7002, 18.9801, 22.2178];

% -----------------------------
% Flatten U with indices and sort
% -----------------------------

% Number of modes (rows * columns)
num_rows = size(U,1);
num_cols = size(U,2);
num_modes = num_rows * num_cols;

% Preallocate the modes array
modes = zeros(num_modes, 3);  % each row will store [l, m, u]

% Fill the array
index = 1;
for l_index = 1:num_rows
    for m_index = 1:num_cols
        u_val = U(l_index, m_index);
        modes(index,:) = [l_index-1, m_index, u_val];  % store (l, m, u)
        index = index + 1;
    end
end

% Sort by u-value ascending
[~, sort_idx] = sort(modes(:,3));
modes = modes(sort_idx,:);

% -----------------------------
% Compute E_lm for allowed modes
% -----------------------------
E_struct = struct();  % to store each mode
allowed_modes = [];

for k = 1:size(modes,1)
    l = modes(k,1);
    m = modes(k,2);
    u_val = modes(k,3);
    
    if u_val > V
        break   % stop when u > V
    end
    
    % compute mode
    E_lm = E_lp_modes(l,m,phi,R,a);  
    % store in struct
    field_name = ['E' num2str(l) num2str(m)];
    E_struct.(field_name) = E_lm;
    
    % keep track of field names for later summation
    allowed_modes = [allowed_modes; {field_name}];
end

% -----------------------------
% Random combination
% -----------------------------
N_modes = length(allowed_modes);
coeffs = rand(1,N_modes);  % random weights
E_total = 0;

for k = 1:N_modes
    E_total = E_total + coeffs(k) * E_struct.(allowed_modes{k});
end

% Normalize
E_total = E_total / max(abs(E_total(:)));

% -----------------------------
% Plot intensity
% -----------------------------
figure;
imagesc(x*1e6,y*1e6,abs(E_total).^2)
axis equal tight
xlabel('µm'); ylabel('µm');
title('MMF random superposition of LP modes');
colorbar;


