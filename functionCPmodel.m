function [P_FIX,P_LO,P_BS,P_UE,P_COD,P_DEC,L_BS,P_BT] = functionCPmodel(valueset)
%% Define parameter values for Value set 1
if valueset == 1
    
    P_FIX = 10;
    
    P_LO = 0.2;
    
    P_BS = 0.4;
    
    P_UE = 0.2;
    
    P_COD = 0.1*10^(-9);
    
    P_DEC = 0.8*10^(-9);
    
    L_BS = 75*10^9;
    
    P_BT = 0.25*10^(-9);
    
% Define parameter values for Value set 2
elseif valueset == 2
    
    P_FIX = 5;
    
    P_LO = 0.1;
    
    P_BS = 0.2;
    
    P_UE = 0.1;
    
    P_COD = 0.01*10^(-9);
    
    P_DEC = 0.08*10^(-9);
    
    L_BS = 750*10^9;
    
    P_BT = 0.025*10^(-9);
    
end
