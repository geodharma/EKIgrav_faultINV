%% EKI run for field data: Mersa Matruh fault, NW Egypt
%% vertical fault, forward model sesar_Eq3 (Eq. 3)
%%
%Laby, D. A., Biswas, A., Heriyanto, M. F., Sungkono, S., Bery, A. A., Ahmed, B. M., Nur Rochman, J. P. G., 
% Purwanto, M. S., & Utama, W. (2026). Interpreting the Gravity Signatures of Faults: 
% Efficient 2D Fault Geometry Estimation Using Regularized Ensemble Kalman Inversion 
% with Sensitivity and Uncertainty Assessment. Submitted.
%%
% Same structure as EKI_invfaultm1.m, only data, search space and
% forward model are changed. Units: x, z, x0 in km; K in mGal
clear, close all, clc;
global x Ne gobs Ny nvar Xmin Xmax lambda
rng(1);   % fixed seed so the run can be reproduced
tic

%% Load data
% col 1 = x (km), col 2 = residual gravity anomaly (mGal)
data=readmatrix("Gravity fault(field_data).xlsx","Sheet","Mesra");
data=data(all(~isnan(data),2),:);
x=data(:,1); gobs=data(:,2);

%% Base level
% residual anomaly goes from negative to positive, base level (constant)
% is solved by least squares for each member (see add_regional.m)
reg_mode='const';

%% Search space [K; z; x0]
Xmin=[ 0; 0;  8];
Xmax=[10; 5; 16];
nvar=length(Xmin); 
sd_kf=0.02;          % data error (mGal)
lambda=1e-1;        % regularization in Kalman gain
Ne=100; Niter=100;  % ensemble size and number of iterations
Ny=length(gobs);

%% Initial ensemble
Xe=initialize(Ne,nvar);   % uniform random between Xmin and Xmax
gcalc=zeros(Ne,Ny);


for i = 1:Ne
    gcalc(i,:) = add_regional(sesar_Eq3(Xe(i,:),x),x,gobs,reg_mode);
    fitness(i,1) = objective_function(gcalc(i,:)');
end

% misfit statistics of the initial ensemble
IQRF = iqr(fitness);
bestFitness = min(fitness);
MedianF = median(fitness);
PDM = [Xe fitness];               % archive of accepted models + misfit
P = [0 bestFitness MedianF IQRF]; % history: iter, best, median, IQR of misfit


%start iterative EKI
for iter = 1:Niter

    gcalcul = zeros(Ne,Ny);
    Xec = EnKF(Ne,sd_kf,Xe,gcalc);   % Kalman update of all members
    
    
    for j = 1 : Ne
        gcalcul (j,:) = add_regional(sesar_Eq3(Xec(j,:),x),x,gobs,reg_mode);
        fitnessN(j,1) = objective_function(gcalcul(j,:)');

        % keep the new member only if the misfit decreases
        if fitnessN(j) < fitness(j)
            fitness(j) = fitnessN(j);
            Xe(j,:) = Xec(j,:);
            gcalc(j,:)  = gcalcul(j,:); 
            PDM = [PDM; Xe(j,:) fitness(j)];
        end
    end

    IQRF = iqr(fitness);
    MedianF = median(fitness);
    [bestFitness, indx] = min(fitness);
    P = [P; iter bestFitness MedianF IQRF];

    fprintf('EKI: Iteration = %i  Error = %.5e \n ', iter, bestFitness)

end
t_run=toc;

%% Posterior summary from final ensemble
MAP = median(Xe,1);   % median as point estimate
unc = iqr(Xe,1);      % interquartile range as uncertainty

bestX = Xe(indx,:);   % member with the lowest misfit
gcal  = gcalc(indx,:);                  % includes base level
rmse  = sqrt(mean((gobs'-gcal).^2));    % mGal
fprintf('\nRMSE best = %.4f mGal | %.1f s\n', rmse, t_run);


figure;
plot(x,gobs,'r*')     % observed
hold on
plot(x,gcal,'b-')     % calculated from bestX (+ base level)
xlabel('Distance (km)'); ylabel('\Deltag (mGal)');


% [best  median  IQR], one row per parameter: K z x0
results = [bestX' MAP' unc']

save('results_Mersa_EKI.mat', 'x','gobs','Xmin','Xmax','reg_mode', ...
     'Ne','Niter','lambda','sd_kf','Xe','gcalc','fitness','PDM','P', ...
     'MAP','unc','bestX','gcal','rmse','results');
