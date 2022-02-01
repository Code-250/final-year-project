function [Mmin,Pmin] = functionMinAntennasPower(SINR,betas,Mstart,Pmax)

%Extract number of BSs/cells and SBSs per cell
[L,~,S] = size(betas);

%Initiate the minimum power required to solve the problem (Inf means
%infeasible)
Pmin = Inf;

%Initiate the number of BS antennas
M = Mstart-1; 

%% Solve the problem until we find a solution that is feasible 
while Pmin > Pmax
    
    %Initiate solution
    Pmin = 0;
    MU = zeros(L,S);
    MU_NEW = ones(L,S);
    STIELTJES = zeros(1,L);
    
    %Increase the number of antennas
    M = M+1;
    
    while (max(max(abs(MU-MU_NEW)./MU_NEW))>1e-6 && Pmin<2*Pmax)
        
        %Here we utilize an algorithm from:
        %Lakshminaryana, S., J. Hoydis, M. Debbah, and M. Assaad. 2010.
        %"Asymptotic analysis of distributed multi-cell beamforming". In: 
        %Proc. IEEE PIMRC. IEEE. 2105-2110.
        
        MU = MU_NEW;
        
        for i = 1:L
            STIELTJES(i) = function_stieltjes(M*reshape(betas(i,:,:),[L,S]),MU,M);
        end
        
        for i = 1:L
            for j=1:S
                MU_NEW(i,j)= SINR/(M*betas(i,i,j)*STIELTJES(i));
            end
        end
        
        %Extract the total average transmit power per BS
        Pmin = 10*log10(sum(sum(MU))/L);
        
    end
    

    
end

Mmin = M;



function m = function_stieltjes(SIGMA,MU,M)

[L,S] = size(MU);
m = 0;
m_new = 1/M;

while abs(m-m_new)/m>1e-6
    
    m = m_new;
    m_new = 1/(1/M*sum(sum((SIGMA.*MU)./(ones(L,S)+SIGMA.*MU*m)))+1);
    
end
