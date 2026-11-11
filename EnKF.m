function Xp=EnKF(Ne,sd_kf,xf,yf)
global gobs Ny Xmin Xmax nvar

Cd = sd_kf^2 * eye(Ny);
mxf=mean(xf,1);
myf=mean(yf,1);
meanxf = repmat(mxf,Ne,1); % Mean of xf
meanyf = repmat(myf,Ne,1);         % Mean of yf
Cxy = (xf-meanxf)'*(yf-meanyf)/(Ne-1);          % Cross-covariance between xf and yf
Cyy = (yf-meanyf)'*(yf-meanyf)/(Ne-1);          % Auto-covariance of yf



lambda = 1e-3;  % regularization parameter
A = Cyy + Cd + lambda * eye(size(Cyy));

K = Cxy / A;  % Calculate the Kalman gain


yobsn=repmat(gobs,1,Ne)+randn(Ny,Ne).*repmat(sd_kf,1,Ne);
X = zeros(Ne, nvar);
for ii=1:Ne
     
    l=rand;

    dX= l*K*(yobsn(:,ii)-yf(ii,:)');             % Calculate the Kalman update
    X(ii,1:nvar)=xf(ii,1:nvar)+dX';  
    Xp(ii,1:nvar)=bound(X(ii,1:nvar),Xmin,Xmax); % Compute the ith proposal
end
