
function E_lm = E_lp_modes(l, m, phi, R, a)

    % Obtaining the corresponding Bessel zero
    u_guess = (m + l/2 - 0.25) * pi;
    u_lm = fzero(@(u) besselj(l, u), u_guess);

    % Operating the complex electric field for each mode
    if l == 0
        E_lm = besselj(0, u_lm * R / a) .* (R <= a);
    else
        E_lm = besselj(l, u_lm * R / a) .* cos(l * phi) .* (R <= a);
    end

end


