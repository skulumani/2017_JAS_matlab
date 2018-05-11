% 26 October 2015
% load the transfer mat files and plot the minimum from each

clc
close all
clear all
addpath(genpath('./'))
constants = crtbp_constants;
constants.control_switch = 'off';

options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);

set(0,'DefaultAxesFontSize',22);

% initialize some plots to use later
traj_fig= figure(1);
hold on
xlabel('$x$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
ylabel('$y$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
title('Trajectory','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

poincare_fig = figure(2);
hold all
grid on
xlabel('$x$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
ylabel('$\dot{x}$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
title('Poincare Section','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

control_fig = figure(3);
hold all
grid on
xlabel('t (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
ylabel('$u$ (N)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
title('Control Input','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

jacobi_fig = figure();
hold all
grid on 
xlabel('$t$ (nondim)', 'interpreter', 'latex', 'FontUnits', 'points', 'FontSize', 22, 'FontName', 'Times')
ylabel('$E$ (nondim)', 'interpreter', 'latex', 'FontUnits', 'points', 'FontSize', 22, 'FontName', 'Times')
title('Jacobi Energy', 'interpreter', 'latex', 'FontUnits', 'points', 'FontSize', 22, 'FontName', 'Times')

% figure for each stage
stage_traj_figs = [];
stage_poincare_figs = [];
% stage_control_figs = [];
for ii = 1:8

    stf = figure();
    hold all
    grid on
    xlabel('$x$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    ylabel('$y$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    title(sprintf('Trajectory Stage %d', ii),'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

    spf = figure();
    hold all
    grid on
    xlabel('$x$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    ylabel('$\dot{x}$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    title(sprintf('Poincare Section Stage %d', ii),'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

    % scf = figure();
    % hold all
    % grid on
    % xlabel('t (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    % ylabel('$u$ (N)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    % title(sprintf('Control Input Stage %d', ii),'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    
    stage_traj_figs = [stage_traj_figs; stf];
    stage_poincare_figs = [stage_poincare_figs; spf];
    % stage_control_figs = [stage_control_figs; scf];
end

%% INITIAL load initial condition (geostationary orbit)


%% TARGET load the stable manifold ( target of the transfer)
manifold_poincare = manifold_parse(traj_fig, poincare_fig);
for ii = 1:8
    temp = manifold_parse(stage_traj_figs(ii), stage_poincare_figs(ii));
end
% generate the target Poincare section
%% loop over the mat files that hold the reachability transfer iterations
reach_switch_vec = {'min_dist' 'min_dist' 'min_dist' 'min_dist' 'min_dist' 'min_dist' 'min_dist' 'max_x'}; 

% concatonate the states into one big trajectory
state_all = [];
control_all = [];
time_all = [];
max_time = 0;
for iter = 1:7
    filename = ['geo_transfer_' num2str(iter) '.mat'];
    load(filename);
    [min_reach, min_man, min_traj,min_costate,min_control,min_time,reach_struct] = minimum_reach(sol_output,manifold_poincare, reach_switch_vec{iter});

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

    text(min_reach(1),min_reach(3),num2str(iter), 'color', 'red', 'fontsize', 12, 'verticalalignment', 'bottom','horizontalalignment','right')
    
    % plot all of the reach states onto their specific plot
    set(0,'CurrentFigure',stage_poincare_figs(iter))
    for jj = 1:length(reach_struct)
        plot(reach_struct(jj).reach_end(1), reach_struct(jj).reach_end(3), 'rs', 'markersize', 8, 'linewidth', 2)
    end
    
    % plot(manifold_poincare(:, 1), manifold_poincare(:, 3), 'g.', 'markersize', 20)
    plot(min_reach(1),min_reach(3),'k.','Markersize',20);
    
    if iter == 1
        plot_trajectories(min_time,[state_all; min_traj(:, 1:4)], constants.e_desired, stage_traj_figs(iter), constants)
    else
        set(0, 'CurrentFigure', stage_traj_figs(iter))
        plot(state_all(:, 1), state_all(:, 2), 'r')
        plot_trajectories(min_time, min_traj(:, 1:4), constants.e_desired, stage_traj_figs(iter), constants)
    end
    %     if iter == 3
%         set(0,'CurrentFigure',poincare_fig)
%         for ii = 1:length(reach_struct)
%             if reach_struct(ii).reach_end(1) > 0.13
%                 plot(reach_struct(ii).reach_end(1),reach_struct(ii).reach_end(3),'r.','markersize',20)
%                 line([min_reach(1) min_man(1)],[min_reach(3) min_man(3)],'linewidth',3)
%             end
%         end
%     end
    
    
    time_all = [time_all;min_time+repmat(max_time,length(min_time),1)];
    max_time = time_all(end);
    state_all = [state_all;min_traj(:,1:4)];
    control_all = [control_all;min_control];
end

% add on the last optimal transfer trajectory to get onto the stable
% manifold
load('geo_transfer_final.mat')
state_all = [state_all;state(:,1:4)];
control_all = [control_all;u];
control_all = control_all * constants.a_scale * constants.sc_mass * constants.km2meter;
time_all = [time_all;t+repmat(max_time,length(t),1)];

% run and print some stats for the final transfer computation
xc0 = [0.17541888429434552000; 0.00000000000000000003; -0.40753131851731217000; 1.71883588181152410000 ]; % geo final transfer working
tspan = 0.355837553671024; % geo transfer final working
xcf = [ 0.175319307882103 -0.000000000000020 -0.282163264918425 2.717676740320596]; % geo transfer final working
[~, ~, ~] = pcrtbp_fixed_tf(xc0, xcf, tspan);

set(0,'CurrentFigure',traj_fig);
plot(state(:,1),state(:,2),'r')

% draw the poincare section black line
line([-constants.mu 0.835],[0 0],'color','k','linewidth',3)

set(0,'CurrentFigure',control_fig);
plot(time_all,control_all)
control_legend = legend('$u_x$','$u_y$');
set(control_legend,'interpreter','Latex','FontUnits','points','FontSize',22,'FontName','Times');

set(0,'CurrentFigure',poincare_fig)
plot(state(end,1),state(end,3),'.','Markersize',20);
text(state(end,1),state(end,3),'Final' ,'color', 'black', 'fontsize', 12, 'verticalalignment', 'middle','horizontalalignment','left')

set(0, 'CurrentFigure', stage_poincare_figs(8))
plot(state(end, 1), state(end, 3), 'b*', 'Markersize', 8, 'linewidth', 2);
plot(manifold_poincare(:, 1), manifold_poincare(:, 3), 'g.', 'markersize', 20)
% propogate the last state forward in time until reaching the periodic
% orbit
stable_man = [ 0.175319307882103 -0.000000000000020 -0.282163264918425 2.717676740320596];
opt_trans = state_all(end,1:4);
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 8],opt_trans,options_cross) ;

% set(0,'CurrentFigure',traj_fig);
plot_trajectories(t, state, constants.e_desired, traj_fig, constants)

set(0, 'CurrentFigure', stage_traj_figs(8))
plot(state_all(:, 1), state_all(:, 2), 'r')
plot_trajectories(t, state, constants.e_desired, stage_traj_figs(8), constants)

% compute the jacobi energy for the whole transfer and plot
jacobi_energy = energyconst(state_all, constants.mu);
set(0, 'CurrentFigure', jacobi_fig)
plot(time_all, jacobi_energy)

keyboard
for ii = 1:8
    set(0, 'CurrentFigure', stage_traj_figs(ii))
    set(gcf, 'PaperPositionMode', 'auto');
    title(sprintf('Trajectory Stage %d', ii),'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

    % draw the poincare section black line
    line([-constants.mu 0.835],[0 0],'color','k','linewidth',3)
    % save the figures
    % savefig(sprintf('./journal_figures/geo_transfer/stage%d_trajectory.fig', ii))
    saveas(stage_traj_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_trajectory.eps', ii), 'epsc')
    % export_fig(stage_traj_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_trajectory.eps', ii), '-eps', '-q101')
    % print(stage_traj_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_trajectory.eps', ii), '-depsc', '-r800')

    % now zoom in and save
    set(0, 'CurrentFigure', stage_traj_figs(ii))
    set(gcf, 'PaperPositionMode', 'auto');
    axis([-0.25, 0.25, -0.25, 0.25])
    % savefig(sprintf('./journal_figures/geo_transfer/stage%d_trajectory_zoom.fig', ii))
    saveas(stage_traj_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_trajectory_zoom.eps', ii), 'epsc')
    % export_fig(stage_traj_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_trajectory_zoom.eps', ii), '-eps', '-q101')

    % print(gcf, sprintf('./journal_figures/geo_transfer/stage%d_trajectory_zoom.eps', ii), '-depsc', '-r800')

    set(0, 'CurrentFigure', stage_poincare_figs(ii))
    set(gcf, 'PaperPositionMode', 'auto');
    % savefig(sprintf('./journal_figures/geo_transfer/stage%d_poincare.fig', ii))
    saveas(stage_poincare_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_poincare.eps', ii), 'epsc')
    % export_fig(stage_poincare_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_poincare.eps', ii), '-eps', '-q101')
    % print(stage_poincare_figs(ii), sprintf('./journal_figures/geo_transfer/stage%d_poincare.eps', ii), '-depsc', '-r800')
end
