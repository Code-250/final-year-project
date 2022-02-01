function channelGaindB = functionExampleSetup_backhaul(S,L)

%Set the length in meters of the total square area
squareLength = 1000;

%Number of BSs per dimension
nbrBSsPerDim = sqrt(L);

%Pathloss exponent
alpha = 3.76;

%Average channel gain in dB at a reference distance of 1 meter. Note that
%-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;

%Standard deviation of shadow fading
sigma_sf = 10;

%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

%Compute all nine alternatives of the BS locations when using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);


%Deploy 81 SBSs on a regular grid in each cell
SBSpositions = zeros(81,L);
for b = 1:L
    
    counter = 0;
    
    for i = 1:9
        
        for j = 1:9
            
            counter = counter + 1;
            SBSpositions(counter,b) = interBSDistance*(-0.5 + (i-1)*1/9 + 1i*(-0.5 + (j-1)*1/9))  + BSpositions(b);

        end
        
    end
    
end


%Randomly select K SBSs within each cell
SELECTION = zeros(S,L);
for b=1:L
    tmp = randperm(81);
    SELECTION(:,b) = SBSpositions(tmp(1:S),b);
end

SBSpositions = SELECTION;


%Prepare to store average channel gains (in dB)
channelGaindB = zeros(S,L,L);


%% Go through all the cells
for l = 1:L
    
    
    %Go through all BSs
    for j = 1:L
        
        %Compute the distance from the SBSs in cell l to BS j with a wrap
        %around topology, where the shortest distance between a SBS and the
        %nine different locations of a BS is considered 
        distancesBSj = min(abs( repmat(SBSpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[S 1]) ),[],2);
        
        %Compute average channel gain using the large-scale fading model in
        %(2.3), while neglecting the shadow fading
        channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
        
    end
    
    %Go through all SBSs in cell l and generate shadow fading
    for k = 1:S
        
        %Generate shadow fading realization
        shadowing = sigma_sf*randn(1,1,L);
        channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        
        %Check if another BS has a larger average channel gain to the SBS
        %than BS l
        while channelGainShadowing(l) < max(channelGainShadowing)
            
            %Generate new shadow fading realizations (until all SBSs in 
            %cell l has its largest average channel gain from BS l)
            shadowing = sigma_sf*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            
        end
        
        %Store average channel gains with shadowing fading
        channelGaindB(k,l,:) = channelGainShadowing;
        
    end
    
end

%Permute dimensions to obtain the designated output
channelGaindB = permute(channelGaindB,[3,2,1]);
