% test to compare propogation of the continuous vs discrete state and
% costate equations
% clear all
clc
% close all
% define constants
constants = crtbp_constants;

% choose initial x0, and h0
x0 = rand(1,4);
h0 = rand(1,4);

t0 = 0;
tf = 2; % figure out how to dimensionalize time
step_exp = 5;
num_steps = round(1*10^step_exp);

constants.control_switch = 'off';
constants.um = 1;

% propogate continuous time system using ODE45
ode_st = tic;
[t_ode45,state_ode45] = ode45(@(t,state)pcrtbp_ode_optimal(t,state,constants),linspace(t0,tf,num_steps),[x0 h0],constants.ode_options);
ode_end = toc(ode_st);

% propogate discrete time system using discrete update equations
% vir_st = tic;
% [t_rect, state_rect,H_rect] = pcrtbp_optimal_var(x0,h0,t0,tf, 'rect', num_steps, constants);
% vir_end = toc(vir_st);
vit_st = tic;
[t_trap, state_trap,H_trap] = pcrtbp_optimal_var(x0,h0,t0,tf, 'trap', num_steps, constants);
vit_end = toc(vit_st);
%  costate_bw = costate_update_bw(t0,tf,num_steps,state_rect(:,1:4),state_rect(:,5:8),constants);

fprintf('tf = %5.2f nondim = %5.2e sec = %5.2f yrs\n', tf, tf*constants.t_scale, tf*constants.t_scale/86400/365);
fprintf('N = %5.2e steps  \n', num_steps);


fprintf('ODE run time %5.4f sec\n',ode_end);
% fprintf('VI RECT run time %5.4f sec\n',vir_end);
fprintf('VI TRAP run time %5.4f sec\n',vit_end);

% plot and compare the differences
% compare trajectory components
comp_fig = figure;
% costate_fig = figure;
set(0,'CurrentFigure',comp_fig)
subplot(4,3,1)
title('State')
hold all;grid on
ylabel('x')
plot(t_ode45,state_ode45(:,1),t_trap,state_trap(:,1))
legend('ODE45''VI-TRAP')
subplot(4,3,4)
hold all;grid on
ylabel('y')
plot(t_ode45,state_ode45(:,2),t_trap,state_trap(:,2))
subplot(4,3,7)
hold all;grid on
ylabel('xd')
plot(t_ode45,state_ode45(:,3),t_trap,state_trap(:,3))
subplot(4,3,10)
hold all;grid on
ylabel('yd')
xlabel('t (nondim)')
plot(t_ode45,state_ode45(:,4),t_trap,state_trap(:,4))

set(0,'CurrentFigure',comp_fig)
subplot(4,3,2)
title('Costate')
hold all;grid on
ylabel('hx')
plot(t_ode45,state_ode45(:,5), t_trap,state_trap(:,5))
% plot(t_rect,costate_bw(:,1))
legend('ODE45', 'VI TRAP')
subplot(4,3,5)
hold all;grid on
ylabel('hy')
plot(t_ode45,state_ode45(:,6),t_trap,state_trap(:,6))
% plot(t_rect,costate_bw(:,2))
subplot(4,3,8)
hold all;grid on
ylabel('hxd')
plot(t_ode45,state_ode45(:,7),t_trap,state_trap(:,7))
% plot(t_rect,costate_bw(:,3))
subplot(4,3,11)
hold all;grid on
ylabel('hyd')
xlabel('t (nondim)')
plot(t_ode45,state_ode45(:,8),t_trap,state_trap(:,8))
% plot(t_rect,costate_bw(:,4))

set(0,'CurrentFigure',comp_fig)
subplot(4,3,3)
title('Difference in costates')
hold all;grid on
ylabel('x or hx')
plot(t_ode45,abs(state_ode45(:,5) - state_trap(:,5)))
legend('VI TRAP - ODE45')
subplot(4,3,6)
hold all;grid on
ylabel('y or hy')
plot(t_ode45,abs(state_ode45(:,6) - state_trap(:,6)))
subplot(4,3,9)
hold all;grid on
ylabel('xd or hxd')
plot(t_ode45,abs(state_ode45(:,7) - state_trap(:,7)))
subplot(4,3,12)
hold all;grid on
ylabel('yd or hyd')
xlabel('t (nondim)')
plot(t_ode45,abs(state_ode45(:,8) - state_trap(:,8)))

%
% % compare costate bw and fw
% set(0,'CurrentFigure',costate_fig)
% subplot(4,2,1)
% title('Costate')
% hold all;grid on
% ylabel('hx')
% plot(t_rect,state_rect(:,5))
% plot(t_rect,costate_bw(:,1))
% plot(t_ode45,state_ode45(:,5))
% legend('VI-RECT FW', 'VI RECT BW', 'ODE45')
% subplot(4,2,3)
% hold all;grid on
% ylabel('hy')
% plot(t_rect,state_rect(:,6))
% plot(t_rect,costate_bw(:,2))
% plot(t_ode45,state_ode45(:,6))
% subplot(4,2,5)
% hold all;grid on
% ylabel('hxd')
% plot(t_rect,state_rect(:,7))
% plot(t_rect,costate_bw(:,3))
% plot(t_ode45,state_ode45(:,7))
% subplot(4,2,7)
% hold all;grid on
% ylabel('hyd')
% xlabel('t (nondim)')
% plot(t_rect,state_rect(:,8))
% plot(t_rect,costate_bw(:,4))
% plot(t_ode45,state_ode45(:,8))
% 
% set(0,'CurrentFigure',costate_fig)
% subplot(4,2,2)
% title('Difference between fw and bw')
% hold all;grid on
% ylabel('hx')
% plot(t_rect,abs(state_rect(:,5) - costate_bw(:,1)))
% subplot(4,2,4)
% hold all;grid on
% ylabel('hy')
% plot(t_rect,abs(state_rect(:,6) - costate_bw(:,2)))
% subplot(4,2,6)
% hold all;grid on
% ylabel('hxd')
% plot(t_rect,abs(state_rect(:,7) - costate_bw(:,3)))
% subplot(4,2,8)
% hold all;grid on
% ylabel('hyd')
% xlabel('t (nondim)')
% plot(t_rect,abs(state_rect(:,8) - costate_bw(:,4)))

% H_comp_fig = figure;
% set(0,'CurrentFigure',H_comp_fig)
% subplot(3,1,1)
% title('Hamiltonian')
% hold all;grid on
% ylabel('H_c with continous ')
% plot(t_ode45,H_cont_rect)
%
% subplot(3,1,2)
% hold all;grid on
% ylabel('H_c with discrete ')
% plot(t_rect,H_cont_rect)
% subplot(3,1,3)
% hold all;grid on
% ylabel('H_d with discrete ')
% plot(t_rect, H_rect)
