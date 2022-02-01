function SE = functionPowerOptimization_prodSINR(signal,interference,Pmax,prelogFactor)


%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);


%% Solve geometric program in (7.8) using CVX
cvx_begin gp
cvx_quiet(true); % This suppresses screen output from the solver

variable rho(K,L);
variable c(K,L);

maximize prod(prod(c))

subject to

for j = 1:L
    
    for k = 1:K
        
        if signal(k,j)>0
            %SINR constraints of UE k in cell j
            c(k,j)*(sum(sum(rho.*interference(:,:,k,j))) + 1) <= (rho(k,j)*signal(k,j));
            
            rho(k,j) >= 0;
            
        else
            %This applies if UE k in cell j is inactive
            c(k,j) == 1;
            rho(k,j) >= 0;
            
        end
        
    end
    
    sum(rho(:,j)) <= Pmax;
    
end

cvx_end

%% Analyze the CVX output and prepare the output variables
if isempty(strfind(cvx_status,'Solved')) %The problem was not solved by CVX, for some reason, and we then consider equal power allocation
    rhoSolution = (Pmax/K)*ones(K,L);
else %The problem was solved by CVX
    rhoSolution = rho;
end

%Refine the solution obtained from CVX using the Matlab command fmincon.
%This is particularly important in case CVX fails to solve the problem
A = kron(eye(L),ones(1,K));
B = Pmax*ones(L,1);
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',50000,'MaxIter',5000);
xend = fmincon(@(x) -functionComputesumSE_DL_poweralloc(x,signal,interference,prelogFactor),rhoSolution(:),A,B,[],[],zeros(K*L,1),[],[],options);
rhoBest = reshape(xend,[K L]);

%Compute the SEs using Theorem 4.6
SE = functionComputeSE_DL_poweralloc(rhoBest,signal,interference,prelogFactor);


function sumSE = functionComputesumSE_DL_poweralloc(rho,signal,interference,prelogFactor)

%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);

%Reshape power variable since fmincon optimizes vectors
rho = reshape(rho,[K L]);

%Prepare to save results
SE = zeros(K,L);

for j = 1:L
    
    for k = 1:K
        
        if signal(k,j) > 0 %Check if the UE k in cell j is active
            
            %Compute the SE in Theorem 4.6 using the formulation in (7.1)
            SE(k,j) = prelogFactor*log2((rho(k,j)*signal(k,j)) / (sum(sum(rho.*interference(:,:,k,j))) + 1));
            
        else %If the UE is inactive
            
            SE(k,j) = 0;
            
        end
        
    end
    
end

%Compute the sum SE of all cells
sumSE = sum(SE(:));
