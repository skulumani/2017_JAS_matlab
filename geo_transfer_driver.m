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
constants.optfsolve = optimoptions(@fsolve,'Display','off','TolFun',1e-3,'TolX',1e-3,...
    'MaxIter',5000,'MaxFunEvals',5000, 'Algorithm', 'trust-region-dogleg','Jacobian','on',...
    'DerivativeCheck','on');

constants.control_switch = 'off';
% geo orbit about the Earth
earth_x0 = [1/constants.l_scale*(35786+6378.137) - constants.mu;0;0;3.07*1/constants.v_scale];

% plot over a long timespan and save all the x axis crossings
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 0.4065],earth_x0,options_cross) ;
plot_trajectories(t, state, energyconst(earth_x0',constants.mu), traj_fig, constants)

% generate the stable manifold for the L1 periodic orbit

load './u=05/l1_manifold.mat'
% number of crossings of the x axis to plot
crossing_sel = 3;
% plot manifold poincare section and trajectories
set(0,'CurrentFigure',traj_fig)
hold on

% pick which invariant manifold to plot
for ii = 1:1:constants.manifold_steps
    uspms = L1_us_manifold_pos_state(:,:,ii);
    usnms = L1_us_manifold_neg_state(:,:,ii);
 
    spms = L1_s_manifold_pos_state(:,:,ii);
    snms = L1_s_manifold_neg_state(:,:,ii);
    
    % maximum time along the manifold
    max_t_uspm = L1_manifold.us_manifold_pos_time_cross(crossing_sel,ii);
    [~,uspms_cross_index ] = min(abs(L1_us_manifold_pos_time(:,ii)  - max_t_uspm));
    
    max_t_usnm = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
    [~,usnms_cross_index ] = min(abs(L1_us_manifold_neg_time(:,ii)  - max_t_usnm));
    
    max_t_spm = L1_manifold.s_manifold_pos_time_cross(crossing_sel,ii);
    [~,spms_cross_index ] = min(abs(L1_s_manifold_pos_time(:,ii)  - max_t_spm));
    
    max_t_snm = L1_manifold.s_manifold_neg_time_cross(crossing_sel,ii);
    [~,snms_cross_index ] = min(abs(L1_s_manifold_neg_time(:,ii)  - max_t_snm));
    
    % plot the manifold
%     plot(uspms(1:uspms_cross_index,1),uspms(1:uspms_cross_index,2),'r')
%     plot(usnms(1:usnms_cross_index,1),usnms(1:usnms_cross_index,2),'r')
    
%     plot(spms(1:spms_cross_index,1),spms(1:spms_cross_index,2),'g')
    plot(snms(1:snms_cross_index,1),snms(1:snms_cross_index,2),'g')
end

plot(state_1(:,1),state_1(:,2),'k','linewidth',4)

line([0.8352 1.1],[0 0],'Linewidth',4,'Color','k')
% Plot Poincare Section


% set(0,'CurrentFigure',poincare_fig)
% 
% for ii = 1:5:constants.manifold_steps
%     
%     
%     usnm_cross_time = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
%     usnm_cross_state = L1_manifold.us_manifold_neg_state_cross(crossing_sel,:,ii);
%     
%     plot(usnm_cross_state(1),usnm_cross_state(3),'g.','Markersize',20)
%     text(usnm_cross_state(1),usnm_cross_state(3), num2str(usnm_cross_time))
% end

%% Call the shooting method function to compute the reachability set for this initial condition
[sol_output]= pcrtbp_shooting(earth_x0, 1*cross_t(2));
% load('l1_reach_earth_u05.mat')

%% Analyze the reachable set by plotting it

poincare_reach_compare(sol_output,traj_fig)