function [Df_mat] = linearized_eom_mat(eq_state, mu)

%   Author:
%       - Shankar Kulumani 4 September 2014
%           - corrected error when only inputting a PCRTBP state (4x1)
%           automatically manages both the spatial and planar case
%       - Shankar Kulumani 16 December 2014
%           - error in Uzz (negative sign)
% function that calculates the linearized eom matrix about a inputted
% equlibrium state

%       (rotating coords)
%
%                 L4
%
% -L3------M1--+-----L1--M2--L2-
%
%                 L5

% define some parameters
mu1 = 1-mu; % larger primary
mu2 = mu; % smaller primary

% switch based on the number of input states for both the spatial and
% planar problems

num_states = length(eq_state);

switch num_states
    case 6 % spatial case
        xe = eq_state(1);
        ye = eq_state(2);
        ze = eq_state(3);
        
        % compute distances to the primary bodies
        r1 = sqrt( (xe + mu2)^2 + ye^2 + ze^2); % distance to larger mass m1
        r2 = sqrt( (xe - mu1)^2 + ye^2 + ze^2); % distance to smaller mass m2
        
        % compute the second partial derivatives of the gravitational potential
        Uxx = 1 - mu1 * ( 1/r1^3 - ( 3*((xe + mu2)^2)/r1^5 ) ) - mu2 * (1/r2^3 - 3*(xe - mu1)^2 /r2^5);
        Uyy = 1 - mu1 * ( 1/r1^3 -  3*ye^2/r1^5 ) - mu2 * ( 1/r2^3 -3*ye^2/r2^5);
        Uzz = -mu1 * ( 1/r1^3 - 3*ze^2/r1^5 ) - mu2 * ( 1/r2^3 - 3*ze^2/r2^5 );
        
        Uxy = 3*mu1*ye*(xe+mu2)/r1^5 + 3*mu2*ye*(xe-mu1)/r2^5; Uyx = Uxy;
        Uxz = 3*mu1*ze*(xe+mu2)/r1^5 + 3*mu2*ze*(xe-mu1)/r2^5; Uzx = Uxz;
        Uyz = 3*mu1*ye*ze/r1^5 + 3*mu2*ye*ze/r2^5; Uzy = Uyz;
        
        % combine everything into the large 6x6 matrix
        G = [Uxx Uxy Uxz;...
            Uyx Uyy Uyz;...
            Uzx Uzy Uzz];
        
        S = 2 * [0 1 0;...
            -1 0 0;...
            0  0 0];
        Df_mat = [zeros(3,3) eye(3,3); ...
            G S];
    case 4 % planar case
        xe = eq_state(1);
        ye = eq_state(2);
%         ze = eq_state(3);
        
        % compute distances to the primary bodies
        r1 = sqrt( (xe + mu2)^2 + ye^2 ); % distance to larger mass m1
        r2 = sqrt( (xe - mu1)^2 + ye^2); % distance to smaller mass m2
        
        % compute the second partial derivatives of the gravitational potential
        Uxx = 1 - mu1 * ( 1/r1^3 - ( 3*((xe + mu2)^2)/r1^5 ) ) - mu2 * (1/r2^3 - 3*(xe - mu1)^2 /r2^5);
        Uyy = 1 - mu1 * ( 1/r1^3 -  3*ye^2/r1^5 ) - mu2 * ( 1/r2^3 -3*ye^2/r2^5);
%         Uzz = -mu1 * ( 1/r1^3 + 3*ze^2/r1^5 ) - mu2 * ( 1/r2^3 + 3*ze^2/r2^5 );
        
        Uxy = 3*mu1*ye*(xe+mu2)/r1^5 + 3*mu2*ye*(xe-mu1)/r2^5; Uyx = Uxy;
%         Uxz = 3*mu1*ze*(xe+mu2)/r1^5 + 3*mu2*ze*(xe-mu1)/r2^5; Uzx = Uxz;
%         Uyz = 3*mu1*ye*ze/r1^5 + 3*mu2*ye*ze/r2^5; Uzy = Uyz;
        
        % combine everything into the large 6x6 matrix
        G = [Uxx Uxy;...
            Uyx Uyy];
        
        S = 2 * [0 1 ;...
            -1 0];
        Df_mat = [zeros(2,2) eye(2,2); ...
            G S];
        
end
