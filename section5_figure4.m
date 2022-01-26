%Empty workspace and close figures
close all;
clear;


%% Propagation and hardware parameters

%Communication bandwidth
B = 0.1*10^(6); %1 kHz

%Number of BS antennas
M = 20;

%PA efficiency
mu = 0.4;

%Range of fixed circuit power value per BS (in Watt)
P_FIX = [0 1 10 20];

%Range of SE values
SE = (0:0.000001:20)';

%Select ratio between noise power and beta_0^0
sigma2_beta = 10^(-6*0.1);

%Compute nu_0 in (5.14)
nu_0 = sigma2_beta/mu;


%Prepare to save simulation results
EE = zeros(length(SE),length(P_FIX));
maxEE = zeros(length(P_FIX),1);
maxSE = zeros(length(P_FIX),1);
maxEE_theory = zeros(length(P_FIX),1);
maxSE_theory = zeros(length(P_FIX),1);


%% Go through all CP values
for index1 = 1:length(P_FIX)
    
    
    %Go through range of SE values
    for index2 = 1:length(SE)
        
        %Compute transmit power using (5.12)
        p = (2^SE(index2) - 1)/(M-1)*sigma2_beta;
        
        %Compute EE using (5.11)
        EE(index2,index1) = B*SE(index2)/(p/mu + P_FIX(index1));
        
    end
    
    if P_FIX(index1) == 0
        
        %Use (5.19) to compute EE when P_FIX = 0
        EE(1,index1) = (M-1)*B/(log(2)*nu_0);
        
    end
    
    %Find the EE-maximizing point on each curve
    [max_value, index_max ] = max(EE(:,index1));
    maxEE(index1) = max_value;
    maxSE(index1) = SE(index_max);
    
    %Find the EE-maximizing pair of SE and EE values, using (5.18) and (5.19)
    argument  = (M-1)*P_FIX(index1)/(nu_0*exp(1)) - 1/exp(1);
    maxSE_theory(index1) = (lambertw(argument) + 1)/log(2);
    maxEE_theory(index1) = (M-1)*B*2^(-maxSE_theory(index1))/(nu_0*log(2));
    
end


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(SE,EE(:,1),'k','LineWidth',1);
plot(SE,EE(:,2),'k-.','LineWidth',1);
plot(SE,EE(:,3),'k--','LineWidth',1);
plot(SE,EE(:,4),'k:','LineWidth',1);

set(gca,'YScale','log');
legend('{P_{FIX} = 0 W}','P_{FIX} = 1 W','P_{FIX} = 10 W','P_{FIX} = 20 W','Location','NorthEast');
xlabel('SE [bit/s/Hz]');
ylabel('EE [bit/Joule]');
title('EE and SE tradeoff with P_{FIX} as a variable.')
axis([0 max(SE) 10^2 10^7]);

plot(maxSE,maxEE,'k-o','LineWidth',1);
plot(maxSE_theory,maxEE_theory,'r-o','LineWidth',1,'MarkerFaceColor','r');
