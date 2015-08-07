% July 13 2015
% script to propogate some trajectories to try and find something feasible
% for the AAS conference

% first idea is to look at a transfer from L1 to L2
clear all
clc 
close all

% load constants
constants = crtbp_constants;
constants.optfsolve = optimoptions(@fsolve,'Display','off','TolFun',1e-4,'TolX',1e-4,...
    'MaxIter',5000,'MaxFunEvals',5000, 'Algorithm', 'trust-region-dogleg','Jacobian','on',...
    'DerivativeCheck','on');

[L_points, ~] = libration_points(constants.mu);
constants.center_vec = L_points(1,:);

% % generate a L1 planar orbit
% [x0_L1, T_L1, E_L1, phi_cross_L1] = periodic_orbit_pcrtbp(1, constants.e_desired, constants);
% % generate a L2 planar orbit
% [x0_L2, T_L2, E_L2, phi_cross_L2] = periodic_orbit_pcrtbp(2, constants.e_desired, constants);
% 
% % propogate each planar orbit with no control and plot to visualize it
constants.control_switch = 'off';

h0_i = zeros(4,1); % assuming no control so it doesn't matter
% [t_L1, state_L1,~] = pcrtbp_optimal_var(x0_L1,h0_i,0,2*T_L1, 'trap', num_steps, constants);
% [t_L2, state_L2,~] = pcrtbp_optimal_var(x0_L2,h0_i,0,2*T_L2, 'trap', num_steps, constants);

% plot the planar orbits
traj_fig= figure(1);

% plot_trajectories(t_L1, state_L1(:,1:4), constants.e_desired, fig_handle, constants)
% plot_trajectories(t_L2, state_L2(:,1:4), constants.e_desired, fig_handle, constants)

% moon bounded orbit
% earth_x0 = [1/constants.l_scale*(35786+6378.137) - constants.mu;0;0;3.07*1/constants.v_scale];
moon_x0 = [1.05;0;0;0.35];
% earth_x0 = [0.75;0;0;0.2883];
% find poincare section of moon orbit over many periods
% plot over a long timespan and save all the x axis crossings
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 20],moon_x0,options_cross) ;
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig, constants)
% load L1 reachable set
% load ./u=01/l1_reach_second
% load ./u=01/l1_reach_first.mat
% load ./u=01/l1_reach_third.mat
% load l1_reach_025.mat
% load l1_reach_second_025.mat
% load ./u=05/l1_reach_first.mat

load ./u=05/l1_manifold.mat

% plot the initial periodic orbit
plot(state_1(:,1),state_1(:,2),'k','linewidth',4)
