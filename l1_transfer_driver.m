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

% now plot the final 'optimal' solution with the transfer to the final orbit

% load the final solution
load("./u=05/l1_reach_final.mat")
min_reach = sol_output(1);
min_reach.constants = sol_output(1).constants;
% combine all the control arrays
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

% plot
control_fig = figure();
plot(time, control)

% plot the optimal trajectory on the trajectory plot
set(0, 'CurrentFigure', traj_fig);
plot(state(:, 1), state(:, 2), 'g')

set(0, 'CurrentFigure', poincare_fig)
plot(state(end, 1), state(end, 3), 'g.')
