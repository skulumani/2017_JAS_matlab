function [E] = energyconst(state,mu)

% write out top matter for this code

%   Author:
%       - Shankar Kulumani 4 September 2014
%           - logic to modify output based on the number of input states.
%           able to handle both planar and spatial cases.
%           - vectorizing code
%
%       - 9 september 2014 - input vector needs to be nx4
% switch based on the number of input states
% 17 october 14 - only outputs the E (energy constant)
num_states = size(state,2);
% seperate out the states
switch num_states
    case 6 % spatial case
        
        pos = state(:,1:3);
        vel = state(:,4:6);
        
        x = pos(:,1);
        y = pos(:,2);
        z = pos(:,3);
        
        xd = vel(:,1);
        yd = vel(:,2);
        zd = vel(:,3);
        
    case 4 % planar case
        
        pos = state(:,1:2);
        vel = state(:,3:4);
        
        x = pos(:,1);
        y = pos(:,2);
        z = zeros(size(x));
        
        xd = vel(:,1);
        yd = vel(:,2);
        zd = zeros(size(x));
        
end
%the distances

d1 = sqrt((x+ mu).^2 + y.^2 + z.^2);
d2 = sqrt( (x+mu-1).^2+ y.^2 + z.^2);

% gravitational potentional
mu1 = 1-mu;
mu2 = mu;

U = mu1./d1 + mu2./d2 + 1/2.*(x.^2 + y.^2); % different value than next two
% U = -1/2*(x^2+y^2) - mu1/d1 -mu2/d2 -1/2*mu1*mu2;
% U = -1/2*(mu1*d1^2 + mu2*d2^2) - mu1/d1 - mu2/d2;

% kinetic energy
% K = 1/2*(vel*vel');

%Compute the energy integral
% E = K + U;

% Jacobi constant
% C = 2*U - vel'*vel  ; % scalar vector product
C = 2.*U - (xd.^2 + yd.^2 + zd.^2);

E = C./-2;