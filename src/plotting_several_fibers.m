
% -----------------------------
% Global parameters
% -----------------------------
lambda = 1.55e-6;        % wavelength [m]
n_core = 1.450;
n_clad = 1.444;

k0 = 2*pi/lambda;

N = 512;                % grid points
L = 50e-6;              % grid half-size
x = linspace(-L, L, N);
y = x;
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
phi = atan2(Y,X);

% Fiber lengths (only phase matters)
L_fiber = 0.1;

% Bessel mode limits
L_max = 10;
M_max = 10;

% -----------------------------
% Precompute Bessel zeros
% -----------------------------
U = zeros(L_max+1, M_max);

for l = 0:L_max
    for m = 1:M_max
        u_guess = (m + l/2 - 0.25)*pi;
        U(l+1,m) = fzero(@(u) besselj(l,u), u_guess);
    end
end

% Flatten and sort modes
modes = [];
for l = 0:L_max
    for m = 1:M_max
        modes = [modes; l, m, U(l+1,m)];
    end
end
modes = sortrows(modes,3);

% -----------------------------
% Define 6 fiber core radii
% -----------------------------
a_list = [5 10 15 20 25 30]*1e-6;   % meters

figure('Color','w','Position',[100 100 1200 700])

for f = 1:length(a_list)

    a = a_list(f);

    % V-number
    NA = sqrt(n_core^2 - n_clad^2);
    V = 2*pi*a/lambda * NA;

    % -----------------------------
    % Compute allowed modes
    % -----------------------------
    E_struct = struct();
    allowed_modes = {};

    for k = 1:size(modes,1)
        l = modes(k,1);
        m = modes(k,2);
        u_val = modes(k,3);

        if u_val > V
            break
        end

        E_lm = E_lp_modes(l, m, phi, R, a);
        fname = sprintf('E%d%d',l,m);
        E_struct.(fname) = E_lm;
        allowed_modes{end+1} = fname;
    end

    % -----------------------------
    % Random superposition
    % -----------------------------
    N_modes = length(allowed_modes);
    coeffs = 2*rand(1,N_modes)-1;
    E_total = 0;

    for k = 1:N_modes
        name = allowed_modes{k};
    
        l_val = str2double(name(2));
        m_val = str2double(name(3));
    
        idx = find(modes(:,1)==l_val & modes(:,2)==m_val, 1);
        u_val = modes(idx,3);
    
        n_eff = n_core * (1 - (u_val/V).^2 / 2);
        beta  = k0 * n_eff;
    
        E_total = E_total + coeffs(k) * ...
                  E_struct.(name) * exp(1i*beta*L_fiber);
    end

    % Normalize
    E_total = E_total / max(abs(E_total(:)));

    % -----------------------------
    % Plot
    % -----------------------------
    subplot(2,3,f)
    imagesc(x*1e6, y*1e6, abs(E_total).^2)
    axis equal tight
    %colormap hot
    colorbar
    xlabel('µm'); ylabel('µm');

    title(sprintf('a = %.0f µm | Modes = %d | V = %.1f',...
          a*1e6, N_modes, V))
end

sgtitle('Multimode Fiber Random LP Mode Superposition')
