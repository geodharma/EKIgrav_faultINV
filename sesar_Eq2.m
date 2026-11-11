function g = sesar_Eq2(m, x)
%Forward model 2D inclined fault, h -> infinity (Eq. 2 R&B 2021).
%   Delta_g(xi) = K * [ pi + atan((xi - x0)/z + cot(alpha)) ]
%
%   m = [K; z; alpha; x0]
    K      = m(1);
    z      = m(2);
    alpha  = m(3);
    x0     = m(4);

    alpha_rad = deg2rad(alpha);
    cot_a     = cot(alpha_rad);

    z_eff = max(z, 1e-6);

    xi_shift = x - x0;
    g = K .* ( pi + atan( xi_shift ./ z_eff + cot_a ) );
end