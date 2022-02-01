function R = functionRlocalscatteringApprox(M,theta,ASDdeg,antennaSpacing)

%INPUT:
%M              = Number of antennas
%theta          = Nominal angle
%ASDdeg         = Angular standard deviation around the nominal angle
%                 (measured in degrees)
%antennaSpacing = (Optional) Spacing between antennas (in wavelengths)
%
%OUTPUT:
%R              = M x M spatial correlation matrix
%

%Set the antenna spacing if not specified by input
if  nargin < 4
    
    %Half a wavelength distance
    antennaSpacing = 1/2;
    
end


%Compute the ASD in radians based on input
ASD = ASDdeg*pi/180; 


%The correlation matrix has a Toeplitz structure, so we only need to
%compute the first row of the matrix
firstRow = zeros(M,1);


%Go through all the columns of the first row
for column = 1:M
    
    %Distance from the first antenna
    distance = column-1;
    
    %Compute the approximated integral as in (2.24)
    firstRow(column) = exp(1i*2*pi*antennaSpacing*sin(theta)*distance)*exp(-ASD^2/2 * ( 2*pi*antennaSpacing*cos(theta)*distance )^2);
    
end

%Compute the spatial correlation matrix by utilizing the Toeplitz structure
R = toeplitz(firstRow);
