function f = objective_function(d_cal)
global gobs
%Objective Function  Data misfit  (Rao & Biswas, 2021)
% d_cal = calculated anomaly, f = misfit (0 = perfect fit)

    d_obs = gobs(:);   
    d_cal = d_cal(:);

    % normalized by |d_obs| + half of anomaly range
    denom = abs(d_obs) + (max(d_obs) - min(d_obs)) / 2;
    f     = mean(((d_obs - d_cal) ./ denom).^2);
  
   
end
