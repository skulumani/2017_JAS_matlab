% 25 Feb 2015 - script to compare performance of variational integrator and
% ODE45 and normal euler integration
% clear all
close all
clc
% define constants
constants = crtbp_constants;
set(0,'DefaultAxesFontSize',22);
% define initial conditions
% [x0, ~, E_i, ~] = periodic_orbit_pcrtbp(constants.l_point, constants.e_desired, constants);
% x0 = x0';
x0 = [0.75;0;0;0.2883]';
t0 = 0;
tf = 200; % figure out how to dimensionalize time
step_exp = 7;
num_steps = round(1*10^step_exp);
width = 5;
height = 2;
traj_fig = figure('PaperPositionMode','auto');
e_fig = figure('PaperPositionMode','auto');
comp_fig = figure('PaperPositionMode','auto');

t_vec = linspace(t0,tf,num_steps);
h = (tf-t0)/num_steps;
% propogate using ODE45
ode_st = tic;
[t_ode45,state_ode45]=ode45(@(t,state)pcrtbp_ode(t,state,constants.mu),[t0 tf], x0,constants.ode_options);
ode_end = toc(ode_st);
% propogate using variational integrator rectangle
vir_st = tic;
[t_rect, state_rect] = pcrtbp_variational(x0,t0,tf, 'rect', num_steps,constants);
vir_end = toc(vir_st);

% use trapezoidal mode
vit_st = tic;
[t_trap, state_trap] = pcrtbp_variational(x0,t0,tf, 'trap', num_steps,constants);
vit_end = toc(vit_st);

% compare trajectory components
set(0,'CurrentFigure',comp_fig)
% subplot(4,2,1)
% title('ODE45 and VI RECT')
% hold all;grid on
% ylabel('x')
% plot(t_ode45,state_ode45(:,1), t_rect,state_rect(:,1))
% legend('ODE45','VI-RECT')
% subplot(4,2,3)
% hold all;grid on
% ylabel('y')
% plot(t_ode45,state_ode45(:,2), t_rect,state_rect(:,2))
% subplot(4,2,5)
% hold all;grid on
% ylabel('xd')
% plot(t_ode45,state_ode45(:,3),t_rect,state_rect(:,3))
% subplot(4,2,7)
% hold all;grid on
% ylabel('yd')
% xlabel('t (nondim)')
% plot(t_ode45,state_ode45(:,4),t_rect,state_rect(:,4))

% plot the trajectories
subplot(4,1,1)
% title('ODE45 and VI TRAP')
hold all;grid on
ylabel('$x$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
plot(t_ode45,state_ode45(:,1),t_trap,state_trap(:,1))
leg=legend('RK45','VI TRAP');
set(leg,'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');
subplot(4,1,2)
hold all;grid on
ylabel('$y$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
plot(t_ode45,state_ode45(:,2),t_trap,state_trap(:,2))
subplot(4,1,3)
hold all;grid on
ylabel('$\dot{x}$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
plot(t_ode45,state_ode45(:,3),t_trap,state_trap(:,3))
subplot(4,1,4)
hold all;grid on
ylabel('$\dot{y}$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
xlabel('$t$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
plot(t_ode45,state_ode45(:,4),t_trap,state_trap(:,4))

% plot the difference
% subplot(4,1,1)
% % title('ODE45 and VI TRAP')
% hold all;grid on
% ylabel('x','interpreter','latex')
% plot(t_ode45,abs(state_ode45(:,1)-state_trap(:,1)))
% legend('RK45','VI TRAP')
% subplot(4,1,2)
% hold all;grid on
% ylabel('y','interpreter','latex')
% plot(t_ode45,abs(state_ode45(:,2)-state_trap(:,2)))
% subplot(4,1,3)
% hold all;grid on
% ylabel('xd','interpreter','latex')
% plot(t_ode45,abs(state_ode45(:,3)-state_trap(:,3)))
% subplot(4,1,4)
% hold all;grid on
% ylabel('yd','interpreter','latex')
% xlabel('t','interpreter','latex')
% plot(t_ode45,abs(state_ode45(:,4)-state_trap(:,4)))

% calculate jacobi energies for both
[E_ode45] = energyconst(state_ode45,constants.mu);
[E_rect] = energyconst(state_rect,constants.mu);
[E_trap] = energyconst(state_trap,constants.mu);

% energy diff
fprintf('tf = %5.2f nondim = %5.2e sec = %5.2f yrs\n', tf, tf*constants.t_scale, tf*constants.t_scale/86400/365);
fprintf('N = %5.2e steps  \n', num_steps);
fprintf('h = %5.2e nondim = %5.2f sec = %5.2f days\n\n',h,h*constants.t_scale, h*constants.t_scale/86400);

fprintf('ODE energy drift %12.10e\n',mean(abs(E_ode45(end)-E_ode45(1))));
fprintf('VI RECT energy drift %12.10e\n',mean(abs(E_rect-E_rect(1))));
fprintf('VI TRAP Energy drift %12.10e\n',mean(abs(E_trap-E_trap(1))));

fprintf('VI RECT mean E var %12.10e\n',var(E_rect));
fprintf('VI TRAP mean E var %12.10e\n',var(E_trap));

fprintf('ODE run time %5.2f sec\n',ode_end);
fprintf('VI RECT run time %5.2f sec\n',vir_end);
fprintf('VI TRAP run time %5.2f sec\n',vit_end);

% plot trajectories and energy differences
set(0,'CurrentFigure',e_fig)

grid on
hold all
xlabel('$t$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
ylabel('$ \Delta E$','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
plot(t_ode45,E_ode45-E_ode45(1))

plot(t_trap,E_trap-E_trap(1))
h=legend('RK45','VI TRAP');
set(h,'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');

plot_trajectories(t_ode45, state_ode45, E_ode45(1), traj_fig, constants)
set(0,'CurrentFigure',traj_fig)
plot(state_trap(:,1),state_trap(:,2),'r')
% plot_trajectories(t_rect, state_rect, E_rect(1), traj_fig, constants)
% plot_trajectories(t_trap, state_trap, E_trap(1), traj_fig, constants)

print(e_fig,'-dpsc2', 'energy.eps')
print(traj_fig,'-dpsc2', 'trajectory.eps')
print(comp_fig,'-dpsc2', 'components.eps')

% print(e_fig,sprintf('./variational_integrator/energy_diff_h1%s%d.eps','e',step_exp))
% print(traj_fig,sprintf('./variational_integrator/trajectory_h1%s%d.eps','e',step_exp))
% print(comp_fig,sprintf('./variational_integrator/state_comp_h1%s%d.eps','e',step_exp))
