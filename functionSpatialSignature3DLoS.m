function a = functionSpatialSignature3DLoS(U,varphi,theta,lambda)



%Define the wave vector in (7.11)
k = -2*pi/lambda * [cos(varphi)*cos(theta); sin(varphi)*cos(theta); sin(theta)];

%Compute the spatial signature as in (7.12)
a = transpose(exp(1i* k'*U));
