function Xp=EnKF(Ne,sd_kf,xf,yf)
% One update step of regularized EKI
% Input : Ne    = ensemble size
%         sd_kf = std of observation perturbation
%         xf    = current ensemble (Ne x nvar)
%         yf    = forward response of xf (Ne x Ny)
% Output: Xp    = updated ensemble, already inside [Xmin, Xmax]
global gobs Ny Xmin Xmax nvar lambda

Cd = sd_kf^2 * eye(Ny);   % data error covariance
mxf=mean(xf,1);
myf=mean(yf,1);
meanxf = repmat(mxf,Ne,1); % Mean of xf
meanyf = repmat(myf,Ne,1);         % Mean of yf
Cxy = (xf-meanxf)'*(yf-meanyf)/(Ne-1);          % Cross-covariance between xf and yf
Cyy = (yf-meanyf)'*(yf-meanyf)/(Ne-1);          % Auto-covariance of yf



if isempty(lambda), lambda = 1e-3; end  % regularization parameter (set in main script)
A = Cyy + Cd + lambda * eye(size(Cyy));   % lambda*I keeps A well conditioned

K = Cxy / A;  % Calculate the Kalman gain


% perturbed observations, one noisy copy of gobs for each member
yobsn=repmat(gobs,1,Ne)+randn(Ny,Ne).*repmat(sd_kf,1,Ne);
X = zeros(Ne, nvar);
for ii=1:Ne
     
    l=rand;   % random step length (0-1)

    dX= l*K*(yobsn(:,ii)-yf(ii,:)');             % Calculate the Kalman update
    X(ii,1:nvar)=xf(ii,1:nvar)+dX';  
    Xp(ii,1:nvar)=bound(X(ii,1:nvar),Xmin,Xmax); % Compute the ith proposal
end
