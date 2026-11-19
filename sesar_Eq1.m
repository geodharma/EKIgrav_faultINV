function g = sesar_Eq1(m, x)
%Forward model 2D inclined fault (Eq. 1 R&B 2021).
%   Delta_g(xi) = K * [ pi + atan((xi - x0)/z + cot(alpha))
%                          - atan((xi - x0)/h + cot(alpha)) ]
%
%   m = [K; z; h; alpha; x0]
%   K in mGal, z h x0 in the unit of x, alpha in degree
    K      = m(1);
    z      = m(2);
    h      = m(3);
    alpha  = m(4);
    x0     = m(5);

    alpha_rad = deg2rad(alpha);
    cot_a     = cot(alpha_rad);

    z_eff = max(z, 1e-6);
    h_eff = max(h, 1e-6);

    xi_shift = x - x0;
    term1 = pi;
    term2 = atan( xi_shift ./ z_eff + cot_a );
    term3 = atan( xi_shift ./ h_eff + cot_a );

    g = K .* ( term1 + term2 - term3 );
end
