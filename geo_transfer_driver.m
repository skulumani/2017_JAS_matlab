% 14 October 2015 - Earth GEO to L1 Stable Manifold transfer

clear all
clc 
close all

% initialize figures for plotting
% trajectory
traj_fig= figure(1);

%% determine Earth geo initial condition

% load constants
constants = crtbp_constants;
constants.optfsolve = optimoptions(@fsolve,'Display','off','TolFun',1e-4,'TolX',1e-4,...
    'MaxIter',5000,'MaxFunEvals',5000, 'Algorithm', 'trust-region-dogleg','Jacobian','on',...
    'DerivativeCheck','on');

constants.control_switch = 'off';
% geo orbit about the Earth
earth_x0 = [1/constants.l_scale*(35786+6378.137) - constants.mu;0;0;3.07*1/constants.v_scale];

% plot over a long timespan and save all the x axis crossings
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 1],earth_x0,options_cross) ;
plot_trajectories(t, state, energyconst(earth_x0',constants.mu), traj_fig, constants)

% generate the stable manifold for the L1 periodic orbit

%% Call the shooting method function to compute the reachability set for this initial condition
[sol_output]= pcrtbp_shooting(earth_x0, cross_t(2));
