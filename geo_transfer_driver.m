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


constants.control_switch = 'off';
% geo orbit about the Earth
% earth_x0 = [1/constants.l_scale*(35786+6378.137) - constants.mu;0;0;3.07*1/constants.v_scale];

% from geo_transfer_1 (minimize distance)
% earth_x0 =[-0.151136352926310;-0.000000000000000;0.132040248107148;-2.314529896548967];
% earth_x0 = [-0.15328451061550127000; -0.00000000000000003188; -0.01424737191265374600; -2.30682745789989330000 ]; % minimize x

% from geo_transfer_2 (minimize distance)
% earth_x0 = [ -0.156738396267405;-0.002517691671468;0.300932441892643;-2.188954314249965]; % from geo_transfer_2.mat (go backward) 
% minimize the x 
% earth_x0 = [-0.16285423791487691000; -0.00000000000000000022; 0.10246279989951737000; -2.12678052734467340000];

% from geo_transfer_3 
% earth_x0 = [-0.160964496092083;0.000000000000015;0.368316756495253;-2.108540717951311];
% minimize x
% earth_x0 = [-0.17089738600471038000; 0.00000000000001411844; 0.21416864780528005000; -1.99867784182314860000];

% from geo_transfer_4
% minimize x
% earth_x0 = [-0.17403568524658589000; 0.00000000000000000000; 0.42346166133524876000; -1.92024337294004920000 ];

% from geo_transfer_5
% minimize x
% earth_x0 = [-0.17371914600914912000; 0.00000000000000007427; 0.57514326654477355000; -1.91895779611010320000];

% from geo_transfer_6
% minimize x
earth_x0 = [-0.17320970217602974000; 0.00000000000000002505; 0.65516376266077869000; -1.98023594146026190000 ];

% propogate forward or backward to find the period of the orbit
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
% forward propogation
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 1],earth_x0,options_cross) ;
% backward propogation
% [t,state,cross_t,cross_state,ie] = ode113(@(t,state)bw_pcrtbp_ode(t,state,constants.mu),[-1 0],earth_x0,options_cross) ;

plot_trajectories(t, state, energyconst(earth_x0',constants.mu), traj_fig, constants)

% initial condition for shooting reachability computation
initial_condition = earth_x0;
period = cross_t(3)+0.001;
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
