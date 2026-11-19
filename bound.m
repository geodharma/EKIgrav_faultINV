function X=bound(X,Xmin,Xmax);
% Keep parameters inside [Xmin, Xmax]
% value outside the bound is reflected back, then clipped if still outside
X = X(:)';
Xmin = Xmin(:)';
Xmax = Xmax(:)';

ii_low = X <Xmin;
ii_up = X>Xmax;
X(ii_low)=2*Xmin(ii_low)-X(ii_low);   % reflect at lower bound
X(ii_up)=2*Xmax(ii_up)-X(ii_up);      % reflect at upper bound

X  = min(max(X, Xmin), Xmax);
end
