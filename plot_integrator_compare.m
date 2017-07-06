function plot_integrator_compare(filename)

    load(filename);
    t0 = 0;
    tf = 200; % figure out how to dimensionalize time
    set(0,'DefaultAxesFontSize',22);

    e_fig = figure('PaperPositionMode','auto'); % energy history
    h_fig = figure('PaperPositionMode', 'auto'); % mean \Delta E
    t_fig = figure('PaperPositionMode', 'auto'); % cpu time

    % plot trajectories and energy differences

    % add to plot
    % loop over the changes in the number of steps/step size
    for ii = 1:length(num_steps)
        ns = num_steps(ii);
        t_vec = linspace(t0, tf, ns);
        h = (tf - t0) / ns;
        set(0,'CurrentFigure',e_fig)
        semilogy(t_vec, abs(Ehist_ode45(ii,1:ns)-Ehist_ode45(ii,1)), 'r')
        hold all
        semilogy(t_vec, abs(Ehist_ode4(ii,1:ns) - Ehist_ode4(ii,1)), 'b')
        semilogy(t_vec, abs(Ehist_trap(ii,1:ns)-Ehist_trap(ii,1)), 'g')

        set(0, 'CurrentFigure', h_fig)
        loglog(h, mean(abs(Ehist_ode45(ii,1:ns) - Ehist_ode45(ii,1))), 'ro')
        hold all
        loglog(h, mean(abs(Ehist_ode4(ii,1:ns) - Ehist_ode4(ii,1))), 'bs')
        loglog(h, mean(abs(Ehist_trap(ii,1:ns) - Ehist_trap(ii,1))), 'gx')

        set(0, 'CurrentFigure', t_fig)
        loglog(h, cputime_ode45(ii), 'ro')
        hold all
        loglog(h, cputime_ode4(ii), 'bs')
        loglog(h, cputime_trap(ii), 'gx')

    end
    set(0,'CurrentFigure',e_fig)
    grid on
    xlabel('$t$ (nondim) ','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    ylabel('$ \Delta E$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

    set(0, 'CurrentFigure', h_fig)
    grid on 
    xlabel('$h$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
    ylabel('Mean $ \Delta E$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

    set(0, 'CurrentFigure', t_fig)
    grid on
    xlabel('$h$ (nondim)', 'interpreter', 'latex', 'FontUnits', 'points', 'FontSize', 22, 'FontName', 'Times')
    ylabel('CPU Time (sec)', 'interpreter', 'latex', 'FontUnits', 'points', 'FontSize', 22, 'FontName', 'Times')
    % h=legend('RK45','VI TRAP');
    % set(h,'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');

    % print(e_fig,'-dpsc2', 'energy.eps')
    % print(traj_fig,'-dpsc2', 'trajectory.eps')
    % print(comp_fig,'-dpsc2', 'components.eps')

    % print(e_fig,sprintf('./variational_integrator/energy_diff_h1%s%d.eps','e',step_exp))
    % print(traj_fig,sprintf('./variational_integrator/trajectory_h1%s%d.eps','e',step_exp))
    % print(comp_fig,sprintf('./variational_integrator/state_comp_h1%s%d.eps','e',step_exp))
