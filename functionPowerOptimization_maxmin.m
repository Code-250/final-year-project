function SE = functionPowerOptimization_maxmin(signal,interference,Pmax,prelogFactor)

%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);

%Check which UEs that have non-zero channels, because these ones are
%excluded (they are considered inactive)
nonzero = reshape(signal,[K*L 1]);
nonzero = nonzero(nonzero>0);

%Initalize the gamma-variables in Algorithm 1
rateLower = 0;
rateUpper = log2(1+Pmax*min(nonzero));

%Set the accuracy of the bisection
delta = 0.01;

%Prepare to save the power solution
rhoBest = zeros(K,L);

%Solve the max-min problem by bisection - iterate until the difference
%between the lower and upper points in the interval is smaller than delta
while norm(rateUpper - rateLower) > delta
    
    %Compute the midpoint of the line. Note that we are performing the
    %bisection in the SE domain instead of the SINR domain as in Algorithm
    %1, since we can then specify delta as the SE difference
    rateCandidate = (rateLower+rateUpper)/2; 
    
    %Transform the midpoints into SINR requirements
    gammaCandidate = 2.^(rateCandidate)-1;
    
    %Solve the feasibility problem using CVX
    [feasible,rhoSolution] = functionFeasibilityProblem_cvx(signal,interference,Pmax,K,L,gammaCandidate);
    
    
    %If the problem was feasible, then replace rateLower with
    %gammaCandidate and store rhoSolution as the new best solution.
    if feasible
        rateLower = rateCandidate;
        rhoBest = rhoSolution;
    else
        %If the problem was not feasible, then replace ratePoint with
        %gammaCandidate
        rateUpper = rateCandidate;
    end
    
end

%Compute the SEs using Theorem 4.6
SE = functionComputeSE_DL_poweralloc(rhoBest,signal,interference,prelogFactor);



function [feasible,rhoSolution] = functionFeasibilityProblem_cvx(signal,interference,Pmax,K,L,SINR)
%Solve the linear feasibility problem in Algorithm 1 using CVX, by adding
%an extra variable so that it becomes a minimization problem with better
%properties.

cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver

variable rho(K,L);  %Variable for the K x L power allocation matrix
variable scaling    %Scaling parameter for power constraints

minimize scaling %Minimize the power indirectly by scaling the power constraints

subject to

for j = 1:L
    
    for k = 1:K
        
        if signal(k,j)>0
            
            %SINR constraints
            SINR*(sum(sum(rho.*interference(:,:,k,j))) + 1) - (rho(k,j)*signal(k,j)) <= 0
            
        end
        
        rho(k,j)>=0
        
    end
    
    sum(rho(:,j)) <= scaling*Pmax;
    
end

scaling >= 0; %Power constraints must be positive

cvx_end


%% Analyze the CVX output and prepare the output variables
if isempty(strfind(cvx_status,'Solved')) %Both the power minimization problem and the feasibility problem are infeasible
    feasible = false;
    rhoSolution = [];
elseif scaling>1 %Only the power minimization problem is feasible
    feasible = false;
    rhoSolution = rho;
else %Both the power minimization problem and feasibility problem are feasible
    feasible = true;
    rhoSolution = rho;
end
