% 14 October 2015 - Earth GEO to L1 Stable Manifold transfer

% clear all
clc 
% close all

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

% need to find the x axis crossing
% earth_x0 =[ -0.151249038163141;0.002422084833166;0.083315071939735;-2.313836947425291]; % geo_transfer
% earth_x0 = [-0.151136352926310;-0.000000000000000;0.132040248107148;-2.314529896548967]; % geo_transfer

% earth_x0 = [ -0.156738396267405;-0.002517691671468;0.300932441892643;-2.188954314249965]; % from geo_transfer_2.mat (go backward) 
earth_x0 = [-0.157056176554181;   0.000000000000000;   0.251675376635517;  -2.188792578810977]; % geo_transfer_2

% earth_x0 = [ -0.160482960979190;-0.002583052200096;0.418007347923897;-2.109025047934479]; % from geo_transfer_3 (go backward)
% earth_x0 = [-0.160964496092083;0.000000000000015;0.368316756495253;-2.108540717951311];

% propogate forward or backward to find the period of the orbit
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
% forward propogation
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 0.3],earth_x0,options_cross) ;
% backward propogation
% [t,state,cross_t,cross_state,ie] = ode113(@(t,state)bw_pcrtbp_ode(t,state,constants.mu),[-1 0],earth_x0,options_cross) ;

plot_trajectories(t, state, energyconst(earth_x0',constants.mu), traj_fig, constants)

% initial condition for shooting reachability computation
initial_condition = earth_x0;
period = cross_t(2);
%% generate the stable manifold for the L1 periodic orbit

load './l1_manifold_geo.mat'
% number of crossings of the x axis to plot
crossing_sel = 3;
% plot manifold poincare section and trajectories
set(0,'CurrentFigure',traj_fig)

% store all manifold poincare section states to compare with the
% reachability set
manifold_poincare = zeros(constants.manifold_steps,4);
% plot the stable manifold of L1
for ii = 1:1:constants.manifold_steps
    uspms = L1_us_manifold_pos_state(:,:,ii);
    usnms = L1_us_manifold_neg_state(:,:,ii);
 
    spms = L1_s_manifold_pos_state(:,:,ii);
    snms = L1_s_manifold_neg_state(:,:,ii);
    
    % maximum time along the manifold
%     max_t_uspm = L1_manifold.us_manifold_pos_time_cross(crossing_sel,ii);
%     [~,uspms_cross_index ] = min(abs(L1_us_manifold_pos_time(:,ii)  - max_t_uspm));
%     
%     max_t_usnm = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
%     [~,usnms_cross_index ] = min(abs(L1_us_manifold_neg_time(:,ii)  - max_t_usnm));
%     
%     max_t_spm = L1_manifold.s_manifold_pos_time_cross(crossing_sel,ii);
%     [~,spms_cross_index ] = min(abs(L1_s_manifold_pos_time(:,ii)  - max_t_spm));
%     
%     max_t_snm = L1_manifold.s_manifold_neg_time_cross(crossing_sel,ii);
%     [~,snms_cross_index ] = min(abs(L1_s_manifold_neg_time(:,ii)  - max_t_snm));
    
% find index of crossing state that hits the U1 section y = 0 x < 0
    snm_index = find(L1_manifold.s_manifold_neg_state_cross(:,1,ii) < 0);
    max_t_snm = L1_manifold.s_manifold_neg_time_cross(snm_index(1),ii);
    [~,snms_cross_index ] = min(abs(L1_s_manifold_neg_time(:,ii)  - max_t_snm));
    
    snm_cross_state = L1_manifold.s_manifold_neg_state_cross(snm_index(1),:,ii);
    manifold_poincare(ii,:) = snm_cross_state;
    % plot the manifold
%     plot(uspms(1:uspms_cross_index,1),uspms(1:uspms_cross_index,2),'r')
%     plot(usnms(1:usnms_cross_index,1),usnms(1:usnms_cross_index,2),'r')
    
%     plot(spms(1:spms_cross_index,1),spms(1:spms_cross_index,2),'g')
    set(0,'CurrentFigure',traj_fig)
    plot(snms(1:snms_cross_index,1),snms(1:snms_cross_index,2),'g')
    set(0,'CurrentFigure',poincare_fig1)
    plot(snm_cross_state(1),snm_cross_state(3),'g.','Markersize',20)
end

% plot(state_1(:,1),state_1(:,2),'k','linewidth',4)

% line([0.8352 1.1],[0 0],'Linewidth',4,'Color','k')
% Plot Poincare Section


% 
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
% [sol_output]= pcrtbp_shooting(initial_condition, period);
load('geo_transfer_3.mat')

%% Analyze the reachable set by plotting it

reach_poincare = poincare_reach_compare(sol_output,traj_fig, poincare_fig1);

%% determine the state closest to the target periodic orbit

% loop over each reach state 
dist_mat = zeros(length(reach_poincare),length(manifold_poincare));
for ii = 1:length(reach_poincare)
    for jj = 1:length(manifold_poincare)
   % compute the distance (x vs xdot) to the current manifold state 
        reach_x = reach_poincare(ii,1);
        reach_xd = reach_poincare(ii,3);
        
        man_x = manifold_poincare(jj,1);
        man_xd = manifold_poincare(jj,3);
        
        dist = sqrt((reach_x-man_x)^2 + (reach_xd-man_xd)^2);
        dist_mat(ii,jj) = dist;
    end
end

% find the minimum and the row/column of the minimum distance
[min_row, row_ind] = min(dist_mat);
[min_col, col_ind] = min(min_row);

% pull out the minimum reach state and corresponding minimum manifold
min_reach = reach_poincare(row_ind(col_ind),:);
min_man = manifold_poincare(col_ind,:);

% plot the minimum state
set(0,'CurrentFigure',poincare_fig1)
line([min_reach(1) min_man(1)],[min_reach(3) min_man(3)])

%% plot just the minimum trajectories and the poincare section
sol_output(row_ind(col_ind)).constants = sol_output(1).constants;
reach_poincare = poincare_reach_compare(sol_output(row_ind(col_ind)),min_traj_fig, min_poincare_fig);