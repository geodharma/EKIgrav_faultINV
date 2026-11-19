%% EKI run for field data: Gazelle (Gazal) fault, Aswan, Egypt
%% inclined fault, forward model sesar_Eq1 (Eq. 1)
%%
%Laby, D. A., Biswas, A., Heriyanto, M. F., Sungkono, S., Bery, A. A., Ahmed, B. M., Nur Rochman, J. P. G., 
% Purwanto, M. S., & Utama, W. (2026). Interpreting the Gravity Signatures of Faults: 
% Efficient 2D Fault Geometry Estimation Using Regularized Ensemble Kalman Inversion 
% with Sensitivity and Uncertainty Assessment. Submitted.
%%
% Same structure as EKI_invfaultm1.m, only data, search space and
% forward model are changed. Units: x, z, h, x0 in m; K in mGal; alpha in degree
clear, close all, clc;
global x Ne gobs Ny nvar Xmin Xmax lambda
rng(1);   % fixed seed so the run can be reproduced
tic

%% Load data
% col 1 = x (m), col 2 = gravity anomaly (mGal)
data=readmatrix("Gravity fault(field_data).xlsx","Sheet","Gazelle Fault");
data=data(all(~isnan(data),2),:);
x=data(:,1); gobs=data(:,2);
Lx=max(x)-min(x);

%% Regional
% Gazelle data are absolute anomaly (~48-73 mGal) on a regional gradient,
% regional is solved by least squares for each member (see add_regional.m)
reg_mode='linear';   % 'linear' or 'const'

%% Search space [K; z; h; alpha; x0]
Xmin=[ 0;    1;     1;    1; min(x)-0.1*Lx];
Xmax=[50; 1000; 10000;  179; max(x)+0.1*Lx];   % x0 covers the profile +/-10%
nvar=length(Xmin); 
sd_kf=0.02;          % data error (mGal)
lambda=1e-3;        % regularization in Kalman gain
Ne=200; Niter=1000; % ensemble size and number of iterations
Ny=length(gobs);

%% Initial ensemble
Xe=initialize(Ne,nvar);   % uniform random between Xmin and Xmax
gcalc=zeros(Ne,Ny);


for i = 1:Ne
    gcalc(i,:) = add_regional(sesar_Eq1(Xe(i,:),x),x,gobs,reg_mode);
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
        gcalcul (j,:) = add_regional(sesar_Eq1(Xec(j,:),x),x,gobs,reg_mode);
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

    if iter==1 || mod(iter,50)==0
        fprintf('EKI: Iteration = %i  Error = %.5e \n ', iter, bestFitness)
    end

end
t_run=toc;

%% Remove non-converged members (misfit > 10x median) from posterior
conv = fitness <= 10*median(fitness);
Xe_all = Xe; gcalc_all = gcalc; fitness_all = fitness;
Xe = Xe(conv,:); gcalc = gcalc(conv,:); fitness = fitness(conv);
fprintf('\nnon-converged members removed: %d of %d\n', sum(~conv), numel(conv));

%% Posterior summary from final ensemble
MAP = median(Xe,1);   % median as point estimate
unc = iqr(Xe,1);      % interquartile range as uncertainty

[bestF, indx] = min(fitness);
bestX = Xe(indx,:);   % member with the lowest misfit
gcal  = gcalc(indx,:);                  % includes regional
rmse  = sqrt(mean((gobs'-gcal).^2));    % mGal

% horizontal position of fault segment at depth z and h
x_top = Xe(:,5) - Xe(:,2).*cotd(Xe(:,4));
x_bot = Xe(:,5) - Xe(:,3).*cotd(Xe(:,4));

fprintf('RMSE best = %.4f mGal | Spearman(alpha,x0) = %+.3f | %.1f s\n', ...
        rmse, corr(Xe(:,4),Xe(:,5),'Type','Spearman'), t_run);
fprintf('x_top = %.2f +/- %.2f m | x_bot = %.2f +/- %.2f m\n', ...
        median(x_top), iqr(x_top), median(x_bot), iqr(x_bot));


figure;
plot(x,gobs,'r*')     % observed
hold on
plot(x,gcal,'b-')     % calculated from bestX (+ regional)
xlabel('Distance (m)'); ylabel('\Deltag (mGal)');


% [best  median  IQR], one row per parameter: K z h alpha x0
results = [bestX' MAP' unc']

save(sprintf('results_Gazelle_EKI_%s.mat',reg_mode), 'x','gobs','Xmin','Xmax', ...
     'reg_mode','Ne','Niter','lambda','sd_kf','Xe','gcalc','fitness', ...
     'Xe_all','gcalc_all','fitness_all','conv','PDM','P','MAP','unc', ...
     'bestX','gcal','rmse','x_top','x_bot','results');
