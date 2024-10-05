close all; clearvars; clc;

%% This script tests the minimal samples that are required to discriminate
%% portal hypertension with statistical properties of SHAPE sensitivity
%% obtained from Monte-Carlo simulation.

%% Polydisperse microbubbles:
%%      mean of sensitivity coeffient = 0.013 dB/mmHg
%%      standard deviation of sensitivity coefficient = 0.017 dB/mmHg
%% Monodisperse microbubbles:
%%      mean of sensitivity coeffient = 0.022 dB/mmHg
%%      standard deviation of sensitivity coefficient = 0.008 dB/mmHg

%% Step 1: 1000 portal-hepatic vein pressure difference cases are randomly 
%% generated to meet uniform distribution in the range between 0 and 20 mmHg
%% Step 2: For each case, generate N samples of SHAPE gradient using the
%% statistical property of SHAPE sensitivity coefficient.
%% Step 3: Plot the ROC curves with different N by using mean value of samples
%% as the discriminative parameter to detect portal hypertension.

rng(1);
sens = 0.022; %0.013; %
stdv = 0.008; %0.017; %

M = 1000;
NS = [2,5,8,10]; %[10,30,50,80]; %
AOCLINES = [];
AOCLABLS = [];

for ndx=1:length(NS)
    N = NS(ndx);
    x = rand(1,M) * 20; %mmHg
    s = x>10;
    y = nan(N,M);
    
    for i=1:length(x)
        y(:,i) = (randn(N,1) * stdv + sens) * x(i);
    end
    
    z = mean(y);
    
    [X,Y,T,AUC,OPTROCPT] = perfcurve(s,z,1);
    index = find(X==OPTROCPT(1) & Y==OPTROCPT(2));
    
    disp('AUC,  1-SPES,  SENS,  THREH');
    [AUC, OPTROCPT, T(index)]
    
    labstr = string(['N=',num2str(N),', AUC=',num2str(AUC,'%1.2f'), ...
        ', Sensitivity=', num2str(OPTROCPT(2),'%1.2f'),...
        ', Specificity=', num2str(1-OPTROCPT(1),'%1.2f'),...
        ', Cutoff=', num2str(T(index),'%1.2f')]);
    
    aocline = plot(X,Y,'LineWidth',1);
    
    AOCLINES = [AOCLINES, aocline];
    AOCLABLS = [AOCLABLS, labstr];
    
    grid on;hold on;
    plot(OPTROCPT(1),OPTROCPT(2),'ro');
    xlabel('FPR (1 - Specificity)') ;
    ylabel('TPR (Sensitivity)');
    title('ROC for Portal Hypertension Assessment');
    
end

legend(AOCLINES,AOCLABLS,'Location','southeast');