function [ydot] = energyfcn(x, y, xdot, mu, e_desired)

% write out top matter for this code

%   Author:
%       - Shankar Kulumani 9 September 2014
%           - modified from energy const equation
%       - 12 September - modifying  planar problem and poincare selection
%       logic
%       - 7 Jan 15 - modifying to be more clear for the planar problem
%           - x, y, xdot can be vectors
% function to output the ydot based on the current state and the
% desired energy

%the distances

d1 = sqrt((x+ mu).^2 + y.^2 );
d2 = sqrt( (x+mu-1).^2+ y.^2 );

% gravitational potentional
mu1 = 1-mu;
mu2 = mu;

U = mu1./d1 + mu2./d2 + 1/2.*(x.^2 + y.^2); % different value than next two
% U = -1/2*(x^2+y^2) - mu1/d1 -mu2/d2 -1/2*mu1*mu2;
% U = -1/2*(mu1*d1^2 + mu2*d2^2) - mu1/d1 - mu2/d2;

C = -2*e_desired*ones(size(x));

ydot = sqrt(-C+2.*U - xdot.^2 );


