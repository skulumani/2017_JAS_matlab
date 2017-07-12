% 20170712 - I don't remember how I generated the data for the control plot
% so I need to extract it and then dimensionalize it and then replot

close all
addpath(genpath('../'))
open ./control_input.fig
h = gcf;
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');

constants = crtbp_constants;

lines = dataObjs{2};
times = get(lines, 'XData');
controls = get(lines, 'YData');
tuy_nondim = times{1};
tux_nondim = times{2};
uy_nondime = controls{1};
ux_nondim = controls{2};

% now dimensionalize the control input

sc_mass = 500; % kilogram mass of spacecraft
km2meter = 1000/1; % convert kilometers to meters

ux_dim = ux_nondim * constants.a_scale * km2meter * sc_mass; % control force in Newtons
uy_dim = uy_nondim * constants.a_scale * km2meter * sc_mass; % control force in newtons

set(0,'DefaultAxesFontSize',22);
figure('PaperPositionMode', 'auto');
hold all
grid on
plot(tux_nondim, ux_dim, 'b', 'DisplayName', '$u_x$')
plot(tuy_nondim, uy_dim, 'r', 'DisplayName', '$u_y$')

title('Control Input','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');
xlabel('$t$ (nondim)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');
ylabel('$u$ (N)','interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');
c_legend = legend('show');

set(c_legend,'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times');
