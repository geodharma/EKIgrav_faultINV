function X=initialize(pop,nvar)
% Initial ensemble, uniform random between Xmin and Xmax
% pop  = ensemble size, nvar = number of parameters
% X    = pop x nvar matrix (one member per row)
global Xmax Xmin 
X=zeros(pop,nvar);
for i=1:nvar 
    X(:,i)=Xmin(i)+(Xmax(i)-Xmin(i))*rand(pop,1); 
end 
