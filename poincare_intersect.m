% July 13 2015
% script to propogate some trajectories to try and find something feasible
% for the AAS conference

% first idea is to look at a transfer from L1 to L2
clear all
clc 
close all

% load L1 reachable set
% load ./u=01/l1_reach_second
% load ./u=01/l1_reach_first.mat
% load ./u=01/l1_reach_third.mat
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
traj_fig= figure(1);
hold all
grid on

poincare_fig = figure(2);
hold all
grid on
xlabel('x axis','interpreter','latex')
ylabel('xdot axis','interpreter','latex')
title('Poincare','interpreter','latex')

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
    set(0,'CurrentFigure',traj_fig);
    plot(state(:,1),state(:,2),'r')
    set(0,'CurrentFigure',poincare_fig);
    plot(state(end,1),state(end,3),'r.')
    
end


reach_poincare = cat(1,reach_struct(:).reach_end); % end states after first reach computation from periodic orbit

% plot_trajectories(t_L1, state_L1(:,1:4), constants.e_desired, fig_handle, constants)
% plot_trajectories(t_L2, state_L2(:,1:4), constants.e_desired, fig_handle, constants)

% moon bounded orbit
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
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig, constants)




% plot the x axis crossing for the Moon bounded orbit
set(0,'CurrentFigure',poincare_fig);
plot(cross_state(:,1),cross_state(:,3),'b.')

% figure out which of the reach_struct trajectories are closest to the
% bounded orbits
% array of all the end states

% extract out the target states that are to the left and below
left_cross = cross_state(cross_state(:,1)< (1-constants.mu),:);
ld_cross = left_cross(left_cross(:,3) < 0,:);
[y,index] = unique(ld_cross(:,1));
ld_cross = ld_cross(index,:);

target_x = ld_cross(:,1);
target_xd = ld_cross(:,3);
% find the limits in the min x and max x direction
min_x = target_x(1);
max_x = target_x(end);

grid_target = linspace(min_x,max_x,100);
% use pchip to interpolate the target orbits
p_target = pchip(target_x, target_xd, grid_target);

figure(3)
hold all
plot(target_x, target_xd,'o',grid_target, p_target,'-')
% interpolate the reachable set now
[y,index] = unique(reach_poincare(:,1));
initial_poincare = reach_poincare(index,:);

initial_x = initial_poincare(:,1);
initial_xd = initial_poincare(:,3);
min_x = initial_x(1);
max_x = initial_x(end);
grid_initial = linspace(min_x,max_x,100);
p_initial = pchip(initial_x, initial_xd, grid_initial);

% plot the target poincare section and the interpolated values

figure(3)
hold all
plot(initial_x, initial_xd,'o',grid_initial,p_initial,'-')

% [x_int,xd_int,iout,jout] = intersections(grid_target,p_target,grid_initial,p_initial);
[x_int,xd_int,iout,jout] = intersections(target_x,target_xd,initial_x,initial_xd);
plot(x_int,xd_int,'g*','markersize',10)

% calculate the target state to control towards
[ydot] = energyfcn(x_int, 0, xd_int, constants.mu, energyconst(moon_x0',constants.mu));

min_reach = sol_output(round(jout));
min_reach.constants = sol_output(1).constants;
% pcrtbp_shooting_min(min_reach)

% plot the data from min reache

% plot the control input (costates)
for jj = 1:num_seg
    x_i = min_reach.x_i;
    h_i = min_reach.h_i;
    start_idx = (jj-1)*num_steps/num_seg+1;
    end_idx = start_idx-1+num_steps/num_seg;
    time(start_idx:end_idx) = min_reach.t(jj,:);
    state(start_idx:end_idx,:) = x_i(:,:,jj);
    costate(start_idx:end_idx,:) = h_i(:,:,jj);
    control(start_idx:end_idx,:) = -min_reach.constants.um * costate(start_idx:end_idx,3:4)./repmat(sqrt(sum(costate(start_idx:end_idx,3:4).^2, 2)),1, 2);
end

figure
plot(time, control)
