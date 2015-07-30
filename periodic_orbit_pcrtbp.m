function [x0_out, t_cross_out, E_out, phi_cross_out] = periodic_orbit_pcrtbp(l_point, e_desired, constants)
%PERIODIC_ORBIT_PCRTBP Generate a periodic orbit about L1/L2 with specific energy
%
% [x0_out t_cross_out, E_out, phi_cross_out] = periodic_orbit(mu, l_point, e_desired)
%
% This will generate a periodic orbit in the planar circular restricted
% three body problem (PCRTBP) about either the L1 or L2 lagrange point at a
% user specified Jacobi energy level. The code uses a linear approximation
% to create two initial seed orbits. These seed orbits are used to
% initialize a numerical continuation/bisection iteratiation method to
% iteratively arrive at the desired orbit.
%
%   Inputs:
%       - mu - mass parameter of the CRTBP system. usually defined in terms
%       the mass ratio of the secondary primary (m2/(m1 + m2))
%       - l_point - lagrange point to generate the orbit about (1 or 2)
%       - e_desired - desired energy of periodic orbit. defined as the
%       negative one half of the jacobi constant (E = C/-2)
%
%   Outputs:
%       - x0 - 4x1 initial condition of the periodic orbit in the rotating
%       reference frame
%       - t_cross_out - 1x1 nondimensional half period of the periodic orbit
%       - E_cross_out - 1x1 jacobi energy constant of periodic orbit
%       - phi_cross_out - 4x4 STM at half period of periodic orbit
%
%   Dependencies:
%       - libration_points.m - calculate location of lagrange points
%       - energyconst.m - calculate jacobi energy constants
%       - periodic_orbit_initial.m - generate initial periodic orbit state
%       using linearization
%       - periodic_orbit_diffcorr.m - differential correction to calculate
%       initial state of periodic orbit for the nonlinear system
%
%   Author:
%       - Shankar Kulumani 31 August 2014
%           - initial version of code
%       - Shankar Kulumani 4 September 2014
%           - modifications found while attempting to expand to spatial
%           case. some dependent functions are only for the planar
%           situation. some functions use logic to work on both the
%           planar/spatial case. should work towards implementing a single
%           code base for both
%       - Shankar Kulumani 24 October 2014
%           - modification to energyconst output since I changed that m
%           file
%           - modified the diffcorr_select_pcrtbp variable name
%
% More detailed help is in the <a href="matlab: help
% periodic_orbit_pcrtbp>extended_help">extended help</a>.

mu = constants.mu;
tol = constants.tol;

gam = (mu/3)^(1/3) ; % estimate of distance between secondary and eq pt

Ax_low  = 3.e-2*gam    ; % initial amplitude (1 of 2)
Ax_high  = 2*3e-2* gam          ; % initial amplitude (2 of 2)

% set up a large matrix that will contain all the calculated periodic
% orbits
max_orbits = 50;
per_orbits_initialstates = zeros(max_orbits,4);
per_orbits_period = zeros(max_orbits,1);
per_orbits_energy = zeros(max_orbits,1);
per_orbits_stm = zeros(max_orbits,16);

% generate the first two differentially corrected periodic orbit intial
% states

[initial_state, ~] = periodic_orbit_initial_pcrtbp(mu, l_point, Ax_low);
[x0, cross_time, phi_cross] = periodic_orbit_diffcorr_pcrtbp(initial_state,constants);
E = energyconst(x0',mu);

per_orbits_initialstates(1,:) = x0';
per_orbits_period(1,1) = [cross_time];
per_orbits_energy(1,:) = [E];
per_orbits_stm(1,:) = reshape(phi_cross,1,16);

[initial_state, ~] = periodic_orbit_initial_pcrtbp(mu, l_point, Ax_high);
[x0, cross_time, phi_cross] = periodic_orbit_diffcorr_pcrtbp(initial_state,constants);
E = energyconst(x0',mu);

per_orbits_initialstates(2,:) = x0';
per_orbits_period(2,1) = [cross_time];
per_orbits_energy(2,:) = [E];
per_orbits_stm(2,:) = reshape(phi_cross,1,16);

% generate a family of periodic orbits up to the selected energy level by
% using the two previous orbits as seeds for iteration

% use the two previous periodic orbits to calculate a delta difference
% use numerical continuation and bisection to bracket the desired energy level

iter = 3;
x0_low = per_orbits_initialstates(1,:)';
x0_high = per_orbits_initialstates(2,:)';
E_low = per_orbits_energy(1,:);
E_high = per_orbits_energy(2,:);

energy_diff = 1;

while iter < max_orbits && energy_diff > tol
    bracket_energy = logical((E_low < e_desired && E_high > e_desired) );
    
    fprintf('PerOrbit Iter = %2g \n', iter)
    
    % check if we bracket the desired energy and decide how to move initial
    % state
    % use the desired energy to determine how to adjust the initial state
    if bracket_energy % bracket energy
        delta = (x0_high-x0_low)/2; % form of bisection method here
        initial_state = x0_low + delta;
    elseif E_high < e_desired % use numerical continuation to increase the orbit
        delta = x0_high-x0_low;
        initial_state = x0_high+delta;
    elseif E_low > e_desired % numerical continuation to decrease orbit
        delta = x0_high-x0_low;
        initial_state = x0_low-delta;
    end
    
    [x0, cross_time, phi_cross] = periodic_orbit_diffcorr_pcrtbp(initial_state,constants);
    E = energyconst(x0',mu);
    
    per_orbits_initialstates(iter,:) = x0';
    per_orbits_period(iter,1) = [cross_time];
    per_orbits_energy(iter,:) = [E];
    per_orbits_stm(iter,:) = reshape(phi_cross,1,16);
    
    % decide if new energy is greater or less than desired
    
    if E < e_desired && E < E_low
        if E_low < e_desired
            % bad
            printf('Iterations are diverging! \n')
            keyboard
        else
            x0_low = x0;
            E_low = E;
        end
    elseif E < e_desired && E > E_low
        x0_low = x0;
        E_low = E;
    elseif E > e_desired && E > E_high
        if E_high > e_desired
            printf('Iterations are diverging! \n')
            keyboard % bad place to be
        else
            x0_high = x0;
            E_high = E;
        end
    elseif E > e_desired && E < E_high
        x0_high = x0;
        E_high = E;
    elseif E < e_desired && E > E_low
        x0_low = x0;
        E_low = E;
    elseif E < e_desired && E > E_high
        x0_high = x0;
        E_high = E;
    elseif E > e_desired && E < E_low
        x0_low = x0;
        E_low = E;
    end
    
    % swap high and low if out of place
    if E_low > E_high
        dum_x = x0_low;
        dum_E = E_low;
        
        E_low = E_high;
        x0_low = x0_high;
        
        E_high = dum_E;
        x0_high = dum_x;
    end
    
    energy_diff = abs(E_high-E_low);
    iter = iter+1;
end

% remove zero entries
per_orbits_initialstates(iter:end,:) = [];
per_orbits_period(iter:end,:) = [];
per_orbits_energy(iter:end,:) = [];
per_orbits_stm(iter:end,:) = [];

% output the last entry (end of iteration state)
x0_out = per_orbits_initialstates(end,:)';
t_cross_out = per_orbits_period(end);
E_out = per_orbits_energy(end);
phi_cross_out = reshape(per_orbits_stm(end,:),4,4);

fprintf('Found correct periodic orbit \n')

if constants.diffcorr_plot
    
    % plot the periodic orbit iterations
    figure(1)
    hold on
    grid on
    axis equal
    
    %numerical integration
    options=odeset('RelTol',1e-13,'AbsTol',1e-22);     %set tolerences
    
    for ii = 1: length(per_orbits_period)
        % plot the periodic orbit family
        tspan=[0 per_orbits_period(ii)]; % only need to integrate half of the orbit
%         [t,state]=ode113(@(t,state)pcrtbp_ode(t,state,mu),tspan, per_orbits_initialstates(ii,:),options);
        [~, state] = trajectory_simulate( per_orbits_initialstates(ii,:), tspan,'pcrtbp', constants);
        plot(state(:,1),state(:,2),'r',state(:,1),-state(:,2),'r')
        drawnow
    end
    
end

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

