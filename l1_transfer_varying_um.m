%% vary the maximum control to see effect on reachability set

clear all
close all
clc

um_array = [0.5, 0.6, 0.7, 0.8, 0.9];
color = ['r', 'g', 'b', 'y', 'm'];
marker = ['o', '+', '*', 'x', 's'];

initial_condition = [  0.815614054266804, 0, 0, 0.192227407664904];
reach_time = 1.307478324303006;
reach_time_array = 1.307478324303006 * [0.95, 1.0, 1.05, 1.1];
%% um change
for ii = 1:size(um_array,2)
    sol_output = pcrtbp_shooting(initial_condition, reach_time, um_array(1));
    save(['./data/l1_varying_um/l1_reach_', num2str(um_array(ii)), '.mat'], 'sol_output');
end

%% now loop over different terminal times


for ii = 1:size(reach_time_array, 2)
    sol_output = pcrtbp_shooting(initial_condition, reach_time_array(ii), 0.5);
    save(['./data/l1_varying_tf_um_5/l1_reach_', num2str(reach_time_array(ii)), '.mat'], 'sol_output');
end

fprintf('Done with tf loop\n')

%% plots the reachable set for varying um
close all
clc

traj_fig = figure(1);
hold all;
grid on;

poincare_fig = figure(2);
hold all
grid on
xlabel('x axis', 'interpreter', 'latex')
ylabel('$\dot{x}$', 'interpreter', 'latex')
title('Poincare', 'interpreter', 'latex')
keyboard
% plot all the reachable set on the Poincare section
for ii = 1:size(um_array, 2)
    load(['./data/l1_varying_um/l1_reach_', num2str(um_array(ii)), '.mat']);
    plot_sol_output_um(sol_output, traj_fig, poincare_fig, color(ii), marker(ii))
    clear sol_output
end
constants = sol_output(1).constants;
% generate the bounded moon orbit
moon_x0 = [1.05;0;0;0.35];
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 20],moon_x0,options_cross) ;
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig, constants)

% plot the moon poincare crossing
set(0,'CurrentFigure',poincare_fig);
plot(cross_state(:,1),cross_state(:,3),'b.')

%% plots the reachable set for varying tf
close all
clc

traj_fig = figure(1);
hold all;
grid on;

poincare_fig = figure(2);
hold all
grid on
xlabel('x axis', 'interpreter', 'latex')
ylabel('$\dot{x}$', 'interpreter', 'latex')
title('Poincare', 'interpreter', 'latex')
keyboard
% plot all the reachable set on the Poincare section
for ii = 1:size(reach_time_array, 2)
    load(['./data/l1_varying_tf/l1_reach_', num2str(reach_time_array(ii)), '.mat']);
    plot_sol_output_tf(sol_output, traj_fig, poincare_fig, color(ii), marker(ii))
end
constants = sol_output(1).constants;
% generate the bounded moon orbit
moon_x0 = [1.05;0;0;0.35];
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 20],moon_x0,options_cross) ;
plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig, constants)

% plot the moon poincare crossing
set(0,'CurrentFigure',poincare_fig);
plot(cross_state(:,1),cross_state(:,3),'b.')

function plot_sol_output_tf(sol_output, traj_fig, poincare_fig, color, marker)
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
    plot(state(:,1),state(:,2),'Color', color, 'Marker', marker)
    set(0,'CurrentFigure',poincare_fig);
    plot(state(end,1),state(end,3),'Color', color, 'Marker', marker)
    
end

end

function plot_sol_output_um(sol_output, traj_fig, poincare_fig, color, marker)
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
    plot(state(:,1),state(:,2),'Color', color, 'Marker', marker)
    set(0,'CurrentFigure',poincare_fig);
    plot(state(end,1),state(end,3),'Color', color, 'Marker', marker)
    
end

end