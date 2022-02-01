function [Hhat_MMSE,C_MMSE,tau_p,R,H,Hhat_EW_MMSE,C_EW_MMSE,Hhat_LS,C_LS] = functionChannelEstimates(R,channelGaindB,nbrOfRealizations,M,K,L,p,f)
%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(M,nbrOfRealizations,K,L,L)+1i*randn(M,nbrOfRealizations,K,L,L));

%Prepare a matrix to save the channel gains per UE
betas = zeros(K,L,L);


%Go through all channels and apply the channel gains to the spatial 
%correlation matrices
for j = 1:L
    
    for l = 1:L
        
        for k = 1:K
            
            if channelGaindB(k,j,l)>-Inf
                
                %Extract channel gain in linear scale
                betas(k,j,l) = 10^(channelGaindB(k,j,l)/10);
                
                %Apply channel gain to correlation matrix
                R(:,:,k,j,l) = betas(k,j,l) * R(:,:,k,j,l);
                
                %Apply correlation to the uncorrelated channel realizations
                Rsqrt = sqrtm(R(:,:,k,j,l));
                H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*H(:,:,k,j,l);
                
            else
                
                betas(k,j,l) = 0;
                R(:,:,k,j,l) = 0;
                H(:,:,k,j,l) = 0;
                
            end
            
        end
        
    end
    
end



%% Perform channel estimation


%Length of pilot sequences
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


%Store identity matrix of size M x M
eyeM = eye(M);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L,f) + 1i*randn(M,nbrOfRealizations,K,L,f));



%Prepare for MMSE estimation

%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M,nbrOfRealizations,K,L,L);

%Prepare to store estimation error correlation matrices
C_MMSE = zeros(M,M,K,L,L);


%Prepare for EW-MMSE estimation
if nargout >= 5
    
    %Prepare to store EW-MMSE channel estimates
    Hhat_EW_MMSE = zeros(M,nbrOfRealizations,K,L,L);
    
    %Prepare to store estimation error correlation matrices
    C_EW_MMSE = zeros(M,M,K,L,L);
    
end

%Prepare for LS estimation
if nargout >= 7
    
    %Prepare to store EW-MMSE channel estimates
    Hhat_LS = zeros(M,nbrOfRealizations,K,L,L);
    
    %Prepare to store estimation error correlation matrices
    C_LS = zeros(M,M,K,L,L);
    
end


%% Go through all cells
for j = 1:L
    
    %Go through all f pilot groups
    for g = 1:f
        
        %Extract the cells that belong to pilot group g
        groupMembers = find(g==pilotPattern)';
        
        %Compute processed pilot signal for all UEs that use these pilots, according to (3.5)
        yp = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
        
        %Go through all UEs
        for k = 1:K
            
            %Compute the matrix that is inverted in the MMSE estimator
            PsiInv = (p*tau_p*sum(R(:,:,k,g==pilotPattern,j),4) + eyeM);
            
            %If EW-MMSE estimation is to be computed
            if nargout >= 5
                %Compute a vector with elements that are inverted in the EW-MMSE estimator
                PsiInvDiag = diag(PsiInv);
            end
            
            %Go through the cells in pilot group g
            for l = groupMembers
                
                %Compute MMSE estimate of channel between BS l and UE k in
                %cell j using (3.9) in Theorem 3.1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = sqrt(p)*RPsi*yp(:,:,k);
                
                %Compute corresponding estimation error correlation matrix, using (3.11)
                C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*RPsi*R(:,:,k,l,j);
                
                
                %If EW-MMSE estimation is to be computed
                if nargout >= 5
                    
                    %Compute EW-MMSE estimate of channel between BS l and
                    %UE k in cell j using (3.33)
                    A_EW_MMSE = diag(sqrt(p)*diag(R(:,:,k,l,j)) ./ PsiInvDiag);
                    Hhat_EW_MMSE(:,:,k,l,j) = A_EW_MMSE*yp(:,:,k);
                    
                    %Compute corresponding estimation error correlation
                    %matrix, using the principle from (3.29)
                    productAR = A_EW_MMSE * R(:,:,k,l,j);
                    
                    C_EW_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_EW_MMSE*PsiInv*A_EW_MMSE';
                
                end
                
                
                %If LS estimation is to be computed
                if nargout >= 7
                    
                    %Compute LS estimate of channel between BS l and UE k
                    %in cell j using (3.35) and (3.36)
                    A_LS = 1/(sqrt(p)*tau_p);
                    Hhat_LS(:,:,k,l,j) = A_LS*yp(:,:,k);
                    
                    %Compute corresponding estimation error correlation
                    %matrix, using the principle from (3.29)
                    productAR = A_LS * R(:,:,k,l,j);
                    
                    C_LS(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_LS*PsiInv*A_LS';
                
                end
                
            end
            
        end

    end
    
end
