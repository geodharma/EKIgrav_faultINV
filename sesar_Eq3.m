function g = sesar_Eq3(m, x)
%FORWARD_RB2021_EQ3 Forward model 2D vertical fault (Eq. 3 R&B 2021).
%   Delta_g(xi) = K * [ 1/2 + (1/pi) * atan((xi - x0)/z) ]
%
%   m = [K; z; x0]
    K  = m(1);
    z  = m(2);
    x0 = m(3);
    z_eff    = max(z, 1e-6);
    xi_shift = x - x0;
    g = K .* ( 0.5 + (1/pi) .* atan( xi_shift ./ z_eff ) );
end