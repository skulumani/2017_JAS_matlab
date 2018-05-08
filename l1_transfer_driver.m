% Driver to generate the plots for L1 to Moon transfer

clear all
clc
close all

% load the reachable set
load("./u=05/l1_reach_first.mat")

constants = sol_output(1).constants;

% plots the reachable set
traj_fig = figure(1);
hold all;
grid on;

poincare_fig = figure(2);
hold all
grid on
xlabel('x axis', 'interpreter', 'latex')
ylabel('$\dot{x}$', 'interpreter', 'latex')
title('Poincare', 'interpreter', 'latex')

% parse out the states from the sol_output
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

% generate the bounded moon orbit
moon_x0 = [1.05;0;0;0.35];
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 20],moon_x0,options_cross) ;
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig, constants)

% plot the moon poincare crossing
set(0,'CurrentFigure',poincare_fig);
plot(cross_state(:,1),cross_state(:,3),'b.')

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

figure(3)
hold all
plot(initial_x, initial_xd,'o',grid_initial,p_initial,'-')

% [x_int,xd_int,iout,jout] = intersections(grid_target,p_target,grid_initial,p_initial);
[x_int,xd_int,iout,jout] = intersections(target_x,target_xd,initial_x,initial_xd);
plot(x_int,xd_int,'g*','markersize',10)

% calculate the target state to control towards
[ydot] = energyfcn(x_int, 0, xd_int, constants.mu, energyconst(moon_x0',constants.mu));

min_reach = sol_output(round(jout));
min_reach.constants = crtbp_constants();

% plot the control input (costates)
for jj = 1:num_seg
    x_i = min_reach.x_i;
    h_i = min_reach.h_i;
    start_idx = (jj-1)*num_steps/num_seg+1;
    end_idx = start_idx-1+num_steps/num_seg;
    min_time(start_idx:end_idx) = min_reach.t(jj,:);
    min_state(start_idx:end_idx,:) = x_i(:,:,jj);
    min_costate(start_idx:end_idx,:) = h_i(:,:,jj);
    min_control(start_idx:end_idx,:) = -0.75 * costate(start_idx:end_idx,3:4)./repmat(sqrt(sum(costate(start_idx:end_idx,3:4).^2, 2)),1, 2);
end

% now define a terminal state to get to
xf_desired = min_state(end, :);
xf_desired(4) = ydot; 
