%% EKI run for synthetic data faults
%%
%Laby, D. A., Biswas, A., Heriyanto, M. F., Sungkono, S., Bery, A. A., Ahmed, B. M., Nur Rochman, J. P. G., 
% Purwanto, M. S., & Utama, W. (2026). Interpreting the Gravity Signatures of Faults: 
% Efficient 2D Fault Geometry Estimation Using Regularized Ensemble Kalman Inversion 
% with Sensitivity and Uncertainty Assessment. Submitted.
%%
clear, close all, clc;
global x Ne gobs Ny nvar Xmin Xmax lambda
rng(1);   % fixed seed so the run can be reproduced
tic

%% Load data
% col 1 = x (m), col 2-4 = gravity anomaly (mGal) of Model 1, 2, 3
data=readmatrix("data_fault_noisyNEW.csv");
x=data(:,1); gobs=data(:,2);   % Model 1 -> column 2

%% True Parameter 
% Model 1: [K; z; h; alpha; x0]
% Model 2: [K; z; alpha; x0]
% Model 3: [K; z; x0]
m_true1 = [   150;     15;    9;     25;       180 ];
%m_true2 = [   100;     15;    45;  180 ];
%m_true3 = [   100;     16;    90];

%% Search space and EKI setting
Xmin=0.5*m_true1;   % lower bound (50% of true model)
Xmax=1.5*m_true1;;  % upper bound (150% of true model)
nvar=length(Xmin); 
sd_kf=0.02;         % std of observation perturbation, Cd = sd_kf^2*I
lambda=1e-3;        % regularization in Kalman gain
Ne=200; Niter=300;  % ensemble size and number of iterations
Ny=length(gobs);

%% Initial ensemble
Xe=initialize(Ne,nvar);   % uniform random between Xmin and Xmax
gcalc=zeros(Ne,Ny);


for i = 1:Ne
    gcalc(i,:) = sesar_Eq1(Xe(i,:),x);
    fitness(i,1) = objective_function(gcalc(i,:)');
end

% misfit statistics of the initial ensemble
IQRF = iqr(fitness);
bestFitness = min(fitness);
best = bestFitness;
MedianF = median(fitness);
PDM = [Xe fitness];               % archive of accepted models + misfit
P = [1 bestFitness MedianF IQRF]; % history: iter, best, median, IQR of misfit


%start iterative EKI
for iter = 1:Niter

    gcalcul = zeros(Ne,Ny);
    Xec = EnKF(Ne,sd_kf,Xe,gcalc);   % Kalman update of all members
    
    
    for j = 1 : Ne
        gcalcul (j,:) =  sesar_Eq1(Xec(j,:),x);
        fitnessN(j,1) = objective_function(gcalcul(j,:)');

        % keep the new member only if the misfit decreases
        if fitnessN(j) < fitness(j);
            fitness(j) = fitnessN(j);
            Xe(j,:) = Xec(j,:);
            gcalc(j,:)  = gcalcul(j,:); 
            PDM = [PDM; Xe(j,:) fitness(j)];
        end
    end

    IQRF = iqr(fitness);
    MedianF = median(fitness);
    [bestFitness, indx] = min(fitness);
    best = bestFitness;
    P = [P; iter bestFitness MedianF IQRF];

    fprintf('EKI: Iteration = %i  Error = %.5f \n ', iter, bestFitness)

end

xpos = PDM; CC=[];


%% Posterior summary from final ensemble
MAP = median(Xe,1);   % median as point estimate
unc = iqr(Xe,1);      % interquartile range as uncertainty

bestX = Xe(indx,:);   % member with the lowest misfit

gcal = sesar_Eq1(bestX,x);


figure;
plot(x,gobs,'r*-')    % observed
hold on
plot(x,gcal,'b-')     % calculated from bestX


% [true  best  median  IQR], one row per parameter
results = [m_true1 bestX' MAP' unc']
