% 30 July 2015
% script to generate the invariant manifolds from the L1 initial orbit

% clear all
clc
close all

constants = crtbp_constants;

% differential correction to generate periodic orbit with specific energy
[x0_1, T_1, ~, ~] = periodic_orbit_pcrtbp(1, constants.e_desired, constants);
x0_i = x0_1;
T_i = T_1;
% x0_i = [0.815614054266804 0 0 0.192227407664904]'; % initial condition on periodic orbit
% T_i = 1.407478324303006; % half period of initial periodic orbit
% constants.poincare_section = 1;
constants.poincare_section = 'transfer'; % x axis section
% constants.poincare_section = 'y_axis';

constants.manifold_plot = 'false';

options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol);
[t_1,state_1] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 2*T_i],x0_i,options_cross) ;
plot_trajectories(t_1, state_1, constants.e_desired, figure(1), constants)

%% manifold globalization
[L1_us_manifold_pos_state, L1_us_manifold_neg_state, L1_s_manifold_pos_state, L1_s_manifold_neg_state, ...
    L1_us_manifold_pos_time, L1_us_manifold_neg_time, L1_s_manifold_pos_time, L1_s_manifold_neg_time, L1_manifold] ...
    = manifold_gen_pcrtbp(x0_i,2*T_i,constants.e_desired,constants);
% % globalize the manifolds from L2 periodic orbit
% [L2_us_manifold_pos_state, L2_us_manifold_neg_state, L2_s_manifold_pos_state, L2_s_manifold_neg_state, ...
%     L2_us_manifold_pos_time, L2_us_manifold_neg_time, L2_s_manifold_pos_time, L2_s_manifold_neg_time, L2_manifold] ...
%     = manifold_gen_pcrtbp(x0_2,2*T_2,constants.e_desired,constants);


save('l1_manifold_geo5',...
    'L1_us_manifold_neg_state','L1_us_manifold_pos_state','L1_us_manifold_neg_time','L1_us_manifold_pos_time',...
    'L1_s_manifold_neg_state','L1_s_manifold_neg_time','L1_s_manifold_pos_state','L1_s_manifold_pos_time',...
    'L1_manifold','t_1','state_1','constants','x0_i','T_i')