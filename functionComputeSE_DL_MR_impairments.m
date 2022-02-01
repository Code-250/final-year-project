function [SE_MR,SE_MR_asymptotic] = functionComputeSE_DL_MR_impairments(channelGainOverNoise,tau_c,M,K,L,p,rho,f,kappatUE,kapparBS,kappatBS,kapparUE)
%Compute DL SE for different transmit precoding schemes using Corollary
%6.6.
%
%INPUT:
%channelGainOverNoise = K x L x L matrix containing the average channel
%                       gains normalized by the noise variance of al
%                       channels (in dB). channelGainOverNoise(k,j,l) is 
%                       the channel gain between between UE k in cell j
%                       and BS l
%tau_c                = Length of coherence block
%M                    = Number of antennas per BS
%K                    = Number of UEs per cell
%L                    = Number of BSs and cells
%p                    = Uplink transmit power per UE (same for everyone)
%rho                  = Downlink transmit power per UE (same for everyone)
%f                    = Pilot reuse factor
%kappatUE             = Hardware quality of the UEs' transmitters
%kapparBS             = Hardware quality of the BSs' receivers
%kappatBS             = Hardware quality of the BSs' transmitters
%kapparUE             = Hardware quality of the UEs' receivers
%
%OUTPUT:
%SE_MR            = K x L matrix where element (k,l) is the downlink SE of
%                   UE k in cell l achieved with MR precoding
%SE_MR_asymptotic = K x L matrix where element (k,l) is the asymptotic
%                   downlink SE of UE k in cell l achieved with MR precoding

%Initiate the hardware qualities as perfect if these are not given as input
if nargin<9
    kappatUE = 1;
end

if nargin<10
    kapparBS = 1;
end

if nargin<11
    kappatBS = 1;
end

if nargin<12
    kapparUE = 1;
end


%Transform the (normalized) average channel gains from dB to linear scale
betaValues = 10.^(channelGainOverNoise/10);

%Compute length of pilot sequences
tau_p = f*K;

%Generate pilot pattern
if f == 1
    
    pilotPattern = ones(L,1);
    
elseif f == 2 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]);
    
elseif f == 4 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]);
    
elseif f == 16 %Only works in the running example with its 16 BSs
    
    pilotPattern = (1:L)';
    
end


%Prepare to store signal and interference terms for all UEs
signalTerm = zeros(K,L);
interferenceTerm = zeros(K,L);

%If the asymptotic SE should be computed
if nargout>1
    
    signalTerm_asymptotic = zeros(K,L);
    interferenceTerm_asymptotic = zeros(K,L);
    
end


%% Go through all cells
for j = 1:L
    
    %Compute the parameter in (6.50)
    G = (1 + kappatBS*(M-1))/(M*kappatBS*kapparUE);
    
    
    %Go through all UEs
    for k = 1:K
        
        %Extract cells that use same pilots as cell j
        groupMembers = find(pilotPattern==pilotPattern(j))';
        
        %Compute the parameter in (6.38)
        psi_jk = 1/( p*kappatUE*kapparBS*tau_p*sum(betaValues(k,groupMembers,j)) + p*(1-kappatUE*kapparBS)*sum(sum(betaValues(:,:,j))) + 1 );
        
        %Compute the numerator of (6.47)
        signalTerm(k,j) = rho*p*(betaValues(k,j,j))^2*tau_p*psi_jk*M;
        
        %Compute parameter in (6.49)
        Flijk = (1+p*betaValues(:,:,j)*psi_jk*( 1 - kappatUE*kapparBS + (1-kappatUE)*kappatBS*kapparBS*(M-1) ))/(kappatBS*kapparUE*kappatUE*kapparBS);
        
        %Compute the non-coherent interference that UE k in cell j causes
        %to all the UEs, which is a part of the first term in the
        %denominator of (6.47)
        interferenceTerm = interferenceTerm + rho*betaValues(:,:,j).*Flijk;
        
        %Compute the coherent interference that UE k in cell j causes to
        %all the pilot-sharing UEs, which is a part of the second term in
        %the denominator of (6.47)
        interferenceTerm(k,groupMembers) = interferenceTerm(k,groupMembers) + rho*p*(betaValues(k,groupMembers,j)).^2*tau_p*psi_jk*M*G;
        
        %Compute the last two terms in the denominator of (6.47)
        interferenceTerm(k,j) = interferenceTerm(k,j) - signalTerm(k,j) + 1/(kappatUE*kapparBS*kappatBS*kapparUE);

        
        %If the asymptotic SE in (6.51) should be computed
        if nargout>1
            
            %Compute the numerator of (6.51), multiplied with psi_jk
            signalTerm_asymptotic(k,j) = rho*(betaValues(k,j,j))^2*psi_jk;
            
            %Compute the part of the first term in the denominator of
            %(6.51) that UE k in cell j contributes to
            interferenceTerm_asymptotic = interferenceTerm_asymptotic + rho*(betaValues(:,:,j)).^2*psi_jk*(1-kappatUE)/(kappatUE*kapparUE*tau_p);

            %Compute the part of the second term in the denominator of
            %(6.51) that UE k in cell j contributes to
            interferenceTerm_asymptotic(k,groupMembers) = interferenceTerm_asymptotic(k,groupMembers) +  rho*(betaValues(k,groupMembers,j)).^2*psi_jk/kapparUE;
            
            %Subtract the numerator from the denominator, which is an
            %alternative way of obtaining the third term in the denominator
            interferenceTerm_asymptotic(k,j) = interferenceTerm_asymptotic(k,j) - rho*(betaValues(k,j,j))^2*psi_jk;
            
        end
        
    end
    
end

%Compute the prelog factor assuming only downlink transmission
prelogFactor = (tau_c-tau_p)/tau_c;

%Compute the final SE expression based on Corollary 6.6
SE_MR = prelogFactor * log2(1+signalTerm ./ real(interferenceTerm));

if nargout>1
    
    %Compute the final asymptotic SE expression based on Corollary 6.7
    SE_MR_asymptotic = prelogFactor * log2(1+signalTerm_asymptotic ./ real(interferenceTerm_asymptotic));
    
end
