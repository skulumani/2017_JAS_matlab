% 14 October 2015 - Earth GEO to L1 Stable Manifold transfer

% clear all
clc 
close all

% initialize figures for plotting
% trajectory
traj_fig= figure(1);
hold on

poincare_fig1 = figure(2);
hold all
grid on
xlabel('x','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
ylabel('$\dot{x}$','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
title('Poincare Section','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')

min_traj_fig = figure(3);
hold all
min_poincare_fig = figure(4);
hold all
%% determine Earth geo initial condition

% load constants
constants = crtbp_constants;

geo_transfer_initial;

constants.control_switch = 'off';

% propogate forward or backward to find the period of the orbit
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
% forward propogation
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 1],earth_x0,options_cross) ;
% backward propogation
% [t,state,cross_t,cross_state,ie] = ode113(@(t,state)bw_pcrtbp_ode(t,state,constants.mu),[-1 0],earth_x0,options_cross) ;

plot_trajectories(t, state, energyconst(earth_x0',constants.mu), traj_fig, constants)

% initial condition for shooting reachability computation
initial_condition = earth_x0;
period = cross_t(3);
%% generate the stable manifold for the L1 periodic orbit

manifold_poincare = manifold_parse(traj_fig, poincare_fig1);

%% Call the shooting method function to compute the reachability set for this initial condition
[sol_output]= pcrtbp_shooting(initial_condition, period);
% load('geo_transfer_4.mat')

%% Analyze the reachable set by plotting it

reach_poincare = poincare_reach_compare(sol_output,traj_fig, poincare_fig1);

%% determine the state closest to the target periodic orbit

[min_reach, min_man, min_traj] = minimum_reach(sol_output,manifold_poincare);

% plot the minimum state
set(0,'CurrentFigure',min_poincare_fig)
line([min_reach(1) min_man(1)],[min_reach(3) min_man(3)])

%% plot just the minimum trajectories and the poincare section
set(0,'CurrentFigure',min_traj_fig)
plot(min_traj(:,1),min_traj(:,2),'.')
