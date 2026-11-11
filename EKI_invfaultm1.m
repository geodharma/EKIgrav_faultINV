%% EKI run for synthetic data faults
%%
%Laby, D. A., Biswas, A., Heriyanto, M. F., Sungkono, S., Bery, A. A., Ahmed, B. M., Nur Rochman, J. P. G., 
% Purwanto, M. S., & Utama, W. (2026). Interpreting the Gravity Signatures of Faults: 
% Efficient 2D Fault Geometry Estimation Using Regularized Ensemble Kalman Inversion 
% with Sensitivity and Uncertainty Assessment. Submitted.
%%
clear, close all, clc;
global x Ne gobs Ny nvar Xmin Xmax 
tic
data=readmatrix("data_fault_noisyNEW.csv");
x=data(:,1); gobs=data(:,2);

%% True Parameter 
m_true1 = [   150;     15;    9;     25;       180 ];
%m_true2 = [   100;     15;    45;  180 ];
%m_true3 = [   100;     16;    90];

Xmin=0.5*m_true1;
Xmax=1.5*m_true1;;
nvar=length(Xmin); 
sd_kf=0.02;
Ne=200; Niter=300;
Ny=length(gobs);

Xe=initialize(Ne,nvar);
gcalc=zeros(Ne,Ny);


for i = 1:Ne
    gcalc(i,:) = sesar_Eq1(Xe(i,:),x);
    fitness(i,1) = objective_function(gcalc(i,:)');
end

IQRF = iqr(fitness);
bestFitness = min(fitness);
best = bestFitness;
MedianF = median(fitness);
PDM = [Xe fitness];
P = [1 bestFitness MedianF IQRF];


%start iterative EKI
for iter = 1:Niter

    gcalcul = zeros(Ne,Ny);
    Xec = EnKF(Ne,sd_kf,Xe,gcalc);
    
    
    for j = 1 : Ne
        gcalcul (j,:) =  sesar_Eq1(Xec(j,:),x);
        fitnessN(j,1) = objective_function(gcalcul(j,:)');

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



MAP = median(Xe,1);
unc = iqr(Xe,1);

bestX = Xe(indx,:);

gcal = sesar_Eq1(bestX,x);


figure;
plot(x,gobs,'r*-')
hold on
plot(x,gcal,'b-')


results = [m_true1 bestX' MAP' unc']
