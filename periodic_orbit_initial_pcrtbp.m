function [x0, T] = periodic_orbit_initial_pcrtbp(mu, l_point, Ax)
%PERIODIC_ORBIT_INITIAL Generate a intial state using linearized EOMs
%
% [initial_state, period] = periodic_orbit_initial(mu, l_point, Ax)
%
% This will generate the initial state of a periodic orbit in the planar
% circular restricted three body problem by using a linearized
% approximation about the specificed lagrange point.
%
%   Inputs:
%       - mu - mass parameter of the CRTBP system. usually defined in terms
%       the mass ratio of the secondary primary (m2/(m1 + m2))
%       - l_point - lagrange point to generate the orbit about (1 or 2)
%       - Ax - amplitude of the periodic orbit in the x direction. will be
%       calculated if not input
%
%   Outputs:
%       - x0 - 4x1 initial condition of the periodic orbit in the rotating
%       reference frame
%       - T_out - 1x1 nondimensional period of the periodic orbit
%
%   Dependencies:
%       - libration_points.m - calculate location of lagrange points
%       - pcrtbp_eig.m - calculate the eigenvalues/eigenvectors for the
%       linearized system and determine the appropriate subspaces
%
%   Author:
%       - Shankar Kulumani 31 August 2014
%           - list revisions
%
% More detailed help is in the <a href="matlab: help
% periodic_orbit_initial>extended_help">extended help</a>.

% mu = constants.mu;
if nargin < 3
gam = (mu/3)^(1/3) ; % estimate of distance between secondary and eq pt

Ax  = 3.e-2*gam    ; % initial amplitude (1 of 2)
end

% determine lagrange points
L_points = libration_points(mu);

switch l_point
    case 1
        L = [L_points(1,:) 0];
    case 2
        L = [L_points(2,:) 0];
    case 3 
        L = [L_points(3,:) 0];
    case 4
        L = [L_points(4,:) 0];
    case 5
        L = [L_points(5,:) 0];
end

eq_state = [L zeros(1,3)]; % only need the x,y state

[Es,Eu,Ec,Vs,Vu,Vc]=pcrtbp_eig(eq_state, mu);

xe = eq_state(1);
mu_bar = mu*abs(xe-1+mu)^(-3) + (1-mu)*abs(xe+mu)^(-3);
% nu = sqrt(-1/2*(mu_bar-2-sqrt(9*mu_bar^2-8*mu_bar)));
nu = abs(imag(Ec(1)));
tau = -(nu^2 + 2*mu_bar+1)/2/nu;

x0= [xe-Ax;0; 0; -Ax*nu*tau];

T= 2*pi/nu;
    
end

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
%
% Examples:
% foo(1,2,3)
%
% See also:
% BAR
% SOMECLASS/SOMEMETHOD
%
%   References
%       - Dynamical systems, the three-body problem and space mission
%       design Koon WS, Lo MW, Marsden JE, Ross SD.  World Scientific;
%       2000.
%

error('This is a placeholder function just for helptext');
end

