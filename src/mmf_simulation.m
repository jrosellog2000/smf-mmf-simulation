
lambda = 1.55e-6;       % wavelength in meters
a = 25e-6;              % core radius in meters
n_core = 1.450;         
n_clad = 1.444;         

NA = sqrt(n_core^2 - n_clad^2);   % numerical aperture
V = 2*pi*a/lambda*NA;             % normalized frequency

N = 512;        % number of grid points
L = 50e-6;      % grid size in meters
x = linspace(-L,L,N);
y = x;
[X,Y] = meshgrid(x,y); 
R = sqrt(X.^2 + Y.^2); 
phi = atan2(Y, X); 

% Matrix of Bessel zeros (rows: l, cols: m)
L_max = 10;   % maximum l
M_max = 10;   % maximum m 

U = zeros(L_max+1, M_max);  % rows: l=0..L_max, columns: m=1..M_max

for l = 0:L_max
    for m = 1:M_max
        % initial value
        u_guess = (m + l/2 - 0.25)*pi;
        U(l+1, m) = fzero(@(u) besselj(l,u), u_guess);
    end
end


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
coeffs = 2*rand(1, N_modes) - 1;   % uniform in [-1, 1]
E_total = 0;

for k = 1:N_modes
    E_total = E_total + coeffs(k) * E_struct.(allowed_modes{k});
end

% Normalize
E_total = E_total / max(abs(E_total(:)));

% -----------------------------
% Display allowed modes
% -----------------------------
fprintf('Fiber V-number: %.4f\n', V);
% Get number of allowed modes
N_modes = length(allowed_modes);
fprintf('\nNumber of allowed LP modes: %d\n', N_modes);

% Display which modes
fprintf('\nAllowed LP modes (l, m, u):\n');
fprintf('  l   m      u_lm\n');
fprintf('--------------------\n');
for k = 1:N_modes
    % Extract l, m, u from field name E_lm
    name = allowed_modes{k};           % e.g., 'E01'
    l_val = str2double(name(2));       % second character
    m_val = str2double(name(3));       % third character
    u_val = modes(modes(:,1)==l_val & modes(:,2)==m_val,3);
    fprintf(' %2d  %2d   %7.4f\n', l_val, m_val, u_val);
end

% -----------------------------
% Plot intensity
% -----------------------------
figure;
imagesc(x*1e6,y*1e6,abs(E_total).^2)
axis equal tight
xlabel('µm'); ylabel('µm');
title('MMF random superposition of LP modes');
colorbar;


