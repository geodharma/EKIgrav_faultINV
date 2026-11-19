function g = add_regional(g, x, gobs, reg_mode)
% Add regional/base level to the calculated anomaly
% the regional is not an EKI parameter, it is solved by least squares
% from the residual (gobs - g) for each member (variable projection)
%   reg_mode = 'none'   -> no regional
%   reg_mode = 'const'  -> regional = c
%   reg_mode = 'linear' -> regional = a + b*x
g = g(:)'; x = x(:); r = gobs(:) - g';

switch reg_mode
    case 'none'
        return
    case 'const'
        A = ones(numel(x),1);
    case 'linear'
        A = [ones(numel(x),1) x];
    otherwise
        error('reg_mode must be ''none'', ''const'' or ''linear''.');
end
g = g + (A*(A\r))';
end
