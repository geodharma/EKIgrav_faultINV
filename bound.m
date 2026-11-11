function X=bound(X,Xmin,Xmax);
X = X(:)';
Xmin = Xmin(:)';
Xmax = Xmax(:)';

ii_low = X <Xmin;
ii_up = X>Xmax;
X(ii_low)=2*Xmin(ii_low)-X(ii_low);
X(ii_up)=2*Xmax(ii_up)-X(ii_up);

X  = min(max(X, Xmin), Xmax);
end