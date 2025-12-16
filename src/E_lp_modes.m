
function E_lm = E_lp_modes(l, m, phi, R, a)

    % Matrix of Bessel zeros
    U = [
        2.4048, 5.5201, 8.6537, 11.7915, 14.9309;
        3.8317, 7.0156, 10.1735, 13.3237, 16.4706;
        5.1356, 8.4172, 11.6198, 14.7960, 17.9598;
        6.3802, 9.7610, 13.0152, 16.2235, 19.4094;
        7.5883, 11.0647, 14.3725, 17.6160, 20.8269;
        8.7715, 12.3386, 15.7002, 18.9801, 22.2178];

    % (1-based indexing)
    row = l + 1;  % row for l
    col = m;      % column = m

    % Obtaining the corresponding Bessel zero
    u_lm = U(row, col);

    % Operating the complex electric field for each mode
    if l == 0
        E_lm = besselj(0, u_lm * R / a) .* (R <= a);
    else
        E_lm = besselj(l, u_lm * R / a) .* cos(l * phi) .* (R <= a);
    end

end


