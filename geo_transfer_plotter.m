% 26 October 2015
% load the transfer mat files and plot the minimum from each

clc
close all

constants = crtbp_constants;
constants.control_switch = 'off';

options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);

% initialize some plots to use later
traj_fig= figure(1);
hold on
xlabel('x','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
ylabel('y','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
title('Trajectory','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')

poincare_fig = figure(2);
hold all
grid on
xlabel('x','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
ylabel('$\dot{x}$','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
title('Poincare Section','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')

%% INITIAL load initial condition (geostationary orbit)


%% TARGET load the stable manifold ( target of the transfer)
manifold_poincare = manifold_parse(traj_fig, poincare_fig);

% generate the target Poincare section
%% loop over the mat files that hold the reachability transfer iterations
reach_switch_vec = {'min_dist' 'min_dist' 'min_dist' 'min_dist' 'min_dist' 'min_dist' 'min_dist' 'max_x'}; 

% concatonate the states into one big trajectory
state_all = [];
for iter = 1:7
    % file name
    % filename = ['./u=0.5_mindist/geo_transfer_' num2str(iter) '.mat'];
%         filename = ['./u=0.5_minx/geo_transfer_minx_' num2str(iter) '.mat'];
%     filename = ['./u=0.75/geo_transfer_' num2str(iter) '.mat'];
    filename = ['geo_transfer_' num2str(iter) '.mat'];
    load(filename);
    % logic to change if we use minimum distance or maximize x value
    
    
    % find the minimum from each to the target set (stable manifold)
    [min_reach, min_man, min_traj] = minimum_reach(sol_output,manifold_poincare, reach_switch_vec{iter});

    
    % take the minimum reach and propogate forward or backward to the
    % nearest x < 0 axis crossing
    if min_reach(1) < 0
        if min_reach(2) < 0 % propogate backwards
            % backward propogation
            [~,~,cross_t,cross_state,ie] = ode113(@(t,state)bw_pcrtbp_ode(t,state,constants.mu),[-1 0],min_reach(1:4),options_cross) ;
        elseif min_reach(2) > 0
            % forward propogation
            [~,~,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 1],min_reach(1:4),options_cross) ;
        end
    elseif min_reach(1) > 0
        if min_reach(2) > 0 % propogate backwards
            % backward propogation
            [~,~,cross_t,cross_state,ie] = ode113(@(t,state)bw_pcrtbp_ode(t,state,constants.mu),[-1 0],min_reach(1:4),options_cross) ;
        elseif min_reach(2) < 0
            % forward propogation
            [~,~,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 1],min_reach(1:4),options_cross) ;
        end
    end
    fprintf('Reach State Iteration %d\n',iter)
    fprintf('earth_x0 = [')
    fprintf('%20.20f; ',cross_state(1,:))
    fprintf(']\n\n')
    
    % plot the minimum Poincare, and trajectory for each iteration
    set(0,'CurrentFigure',traj_fig);
    plot(min_traj(:,1),min_traj(:,2),'k')

    set(0,'CurrentFigure',poincare_fig)
    plot(min_reach(1),min_reach(3),'k.','Markersize',20);
%     line([min_reach(1) min_man(1)],[min_reach(3) min_man(3)])
    text(min_reach(1),min_reach(3),num2str(iter), 'color', 'red', 'fontsize', 12, 'verticalalignment', 'bottom','horizontalalignment','right')
    
    state_all = [state_all;min_traj];
%     keyboard
end

% add on the last optimal transfer trajectory to get onto the stable
% manifold
load('geo_transfer_final.mat')
state_all = [state_all;state(:,1:4)];

set(0,'CurrentFigure',traj_fig);
plot(state(:,1),state(:,2),'r')

set(0,'CurrentFigure',poincare_fig)
plot(state(end,1),state(end,3),'.','Markersize',20);
text(state(end,1),state(end,3),'Final' ,'color', 'black', 'fontsize', 12, 'verticalalignment', 'middle','horizontalalignment','left')
    
% propogate the last state forward in time until reaching the periodic
% orbit
stable_man = [ 0.175319307882103 -0.000000000000020 -0.282163264918425 2.717676740320596];
opt_trans = state_all(end,1:4);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 8],opt_trans,options_cross) ;

% set(0,'CurrentFigure',traj_fig);
plot_trajectories(t, state, constants.e_desired, traj_fig, constants)

% draw the poincare section black line
line([-constants.mu 0.835],[0 0],'color','k','linewidth',3)