% 26 October 2015
% load the transfer mat files and plot the minimum from each

clc
% close all

constants = crtbp_constants;
constants.control_switch = 'off';

options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@events_xcross_nostop);

% initialize some plots to use later
traj_fig= figure(1);
hold on

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
%% loop over the mat files that hold the transfer iterations
for iter = 1:4
    % file name
% filename = ['./u=0.5_mindist/geo_transfer_' num2str(iter) '.mat'];
%     filename = ['./u=0.5_minx/geo_transfer_minx_' num2str(iter) '.mat'];
filename = ['./u=0.75/geo_transfer_' num2str(iter) '.mat'];
    load(filename);
    % find the minimum from each to the target set (stable manifold)
    [min_reach, min_man, min_traj] = minimum_reach(sol_output,manifold_poincare);
    % plot the minimum Poincare, and trajectory for each iteration
    set(0,'CurrentFigure',traj_fig);
    plot(min_traj(:,1),min_traj(:,2),'.')
    
    set(0,'CurrentFigure',poincare_fig)
    plot(min_reach(1),min_reach(3),'.','Markersize',20);
    line([min_reach(1) min_man(1)],[min_reach(3) min_man(3)])
    
    % take the minimum reach and propogate forward or backward to the
    % nearest x axis crossing
    if min_reach(2) < 0 % propogate backwards
        % backward propogation
        [t,state,cross_t,cross_state,ie] = ode113(@(t,state)bw_pcrtbp_ode(t,state,constants.mu),[-1 0],min_reach(1:4),options_cross) ;
    elseif min_reach(2) > 0
        % forward propogation
        [t,state,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,constants.mu),[0 1],min_reach(1:4),options_cross) ;
    end
    
    fprintf('Reach State Iteration %d\n',iter)
    fprintf('earth_x0 = [')
    fprintf('%20.20f; ',cross_state(1,:))
    fprintf(']\n\n')
    
%     keyboard
end

