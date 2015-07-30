% plot the Poincare section to compare the reachability set against thhose
% of the invariant manifold 
% compare times of flight

close all
clc
clear all
%% load u=0.5 data and plot the reachable set on a poincare section
load ./u=05/l1_reach_first.mat
% load l1_reach_earth_u05

constants = sol_output(1).constants;

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
traj_fig1= figure;
hold all
grid on
traj_fig2 = figure();

poincare_fig = figure;
hold all
grid on
xlabel('x','interpreter','latex','FontUnits','points','FontSize',19,'FontName','Times')
ylabel('$\dot{x}$','interpreter','latex','FontUnits','points','FontSize',19,'FontName','Times')
title('Poincare Section','interpreter','latex','FontUnits','points','FontSize',19,'FontName','Times')

% parse out the states for each theta angle
num_steps = sol_output(1).constants.num_steps;
num_seg = sol_output(1).constants.num_seg;
num_theta = length(sol_output);
num_states = sol_output(1).constants.num_states;

reach_struct(num_theta) = struct('state',[],'costate',[],'reach_end',[]);
for ii = 1:num_theta % loop over theta angles (poincare directions)
    state = zeros(num_steps,num_states);
    costate = zeros(num_steps,num_states);
    
    % loop over the segments and combine trajectories into a big array
    for jj = 1:num_seg
       x_i = sol_output(ii).x_i;
       h_i = sol_output(ii).h_i;
       start_idx = (jj-1)*num_steps/num_seg+1;
       end_idx = start_idx-1+num_steps/num_seg;
       state(start_idx:end_idx,:) = x_i(:,:,jj);
       costate(start_idx:end_idx,:) = h_i(:,:,jj);
    end
    
    reach_struct(ii).state= state;
    reach_struct(ii).costate = costate;
    reach_struct(ii).reach_end = [state(end,:) costate(end,:)];
    % plot the trajectories on the same plot of the moon periodic orbit
    set(0,'CurrentFigure',traj_fig1);
    plot(state(:,1),state(:,2),'r')
    set(0,'CurrentFigure',poincare_fig);
    plot(state(end,1),state(end,3),'r.')
    
end


reach_poincare = cat(1,reach_struct(:).reach_end);
%% plot the moon target bounded orbit states
moon_x0 = [1.05;0;0;0.35];
% moon_x0 = [1.02;0;0;0.75];
% earth_x0 = [1/constants.l_scale*(35786+6378.137) - constants.mu;0;0;3.07*1/constants.v_scale];

% find the max from the last iteration
% [v,i]= min(reach_poincare(:,1));
% earth_x0 = sol_output(i).x_i(end,:,end)';

% find poincare section of moon orbit over many periods
% plot over a long timespan and save all the x axis crossings
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 20],moon_x0,options_cross) ;
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig1, constants)

line([0.8352 1.1],[0 0],'Linewidth',4,'Color','k')

% plot the x axis crossing for the Moon bounded orbit
set(0,'CurrentFigure',poincare_fig);
plot(cross_state(:,1),cross_state(:,3),'b.')
%% load the manifold from the periodic orbit
load './u=05/l1_manifold.mat'
crossing_sel = 3;
% plot manifold poincare section and trajectories
set(0,'CurrentFigure',traj_fig2)
hold on

for ii = 1:5:constants.manifold_steps
    uspms = L1_us_manifold_pos_state(:,:,ii);
    usnms = L1_us_manifold_neg_state(:,:,ii);
 
    spms = L1_s_manifold_pos_state(:,:,ii);
    snms = L1_s_manifold_neg_state(:,:,ii);
    
    max_t_usnm = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
    [~,usnms_cross_index ] = min(abs(L1_us_manifold_neg_time(:,ii)  - max_t_usnm));
    
    plot(usnms(1:usnms_cross_index,1),usnms(1:usnms_cross_index,2),'r')
end

plot(state_1(:,1),state_1(:,2),'k','linewidth',4)
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig2, constants)
line([0.8352 1.1],[0 0],'Linewidth',4,'Color','k')
% Plot Poincare Section


set(0,'CurrentFigure',poincare_fig)

for ii = 1:5:constants.manifold_steps
    
    
    usnm_cross_time = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
    usnm_cross_state = L1_manifold.us_manifold_neg_state_cross(crossing_sel,:,ii);
    
    plot(usnm_cross_state(1),usnm_cross_state(3),'g.','Markersize',20)
    text(usnm_cross_state(1),usnm_cross_state(3), num2str(usnm_cross_time))
end

