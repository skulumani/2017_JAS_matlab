% function [C_grid, E_grid] = energy_contour(mu,fig_handle)
% clear all
close all
% clc
constants = crtbp_constants;
fig_handle = figure;
mu = 0.1;
set(0,'DefaultAxesFontSize',22);
switch constants.plot_scale
    case 'dim'
        l_scale = constants.l_scale;
        units_label = '(km)';
    case 'nondim'
        l_scale = 1;
        units_label = '(nondim)';
end

%Contour plot of the zero velocity energy curve
[X,Y]=meshgrid(-1.5:0.01:1.5);
r1=((X + mu).^2+Y.^2 ).^(1/2);
r2=((X + mu-1).^2+Y.^2 ).^(1/2);
C_grid=2*((1/2)*(X.^2 + Y.^2)+(1-mu)./r1+mu./r2); % jacobi constant
E_grid = C_grid/-2;

% Potential Grid
U_grid = 1/2*(X.^2 + Y.^2) + (1-mu)./r1 + mu./r2;

% calculate lagrange points for this mu
[L_points, ~] = libration_points(mu);

L1 = [L_points(1,:) 0];
L2 = [L_points(2,:) 0];
L3 = [L_points(3,:) 0];
L4 = [L_points(4,:) 0];
L5 = [L_points(5,:) 0];

E1 = energyconst([L_points(1,:) 0 0],mu);
E5 = energyconst([L_points(5,:) 0 0],mu);

E_lev = linspace(E1-0.005,E5,10);

primary_1 = [-mu 0 0];
primary_2 = [1-mu 0 0];

% plot contours about the lagrange points
set(0,'CurrentFigure',fig_handle)
hold all
grid on
title_string = ['Energy Surface Contours ($\mu = ' num2str(mu) ')$'];
title(title_string,'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
xlabel(['x' units_label],'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
ylabel(['y' units_label],'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')
zlabel(['z' units_label],'interpreter','latex','FontUnits','points','FontSize',22,'FontName','Times')

contour(X,Y,E_grid,E_lev,'r');

% plot libration points
plot3(l_scale*L1(1), 0,0, 'k.', l_scale*L2(1), 0,0, 'k.', l_scale*L3(1), 0,0,'k.','Markersize',10)
plot3(l_scale*L4(1), l_scale*L4(2),0,'k.','Markersize',10)
plot3(l_scale*L5(1), l_scale*L5(2),0,'k.','Markersize',10)

text(l_scale*L1(1),0, '$L_1$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',12,'FontName','Times')
text(l_scale*L2(1),0, '$L_2$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',12,'FontName','Times')
text(l_scale*L3(1),0, '$L_3$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',12,'FontName','Times')
text(l_scale*L4(1),l_scale*L4(2), '$L_4$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',12,'FontName','Times')
text(l_scale*L5(1),l_scale*L5(2), 'L_5','HorizontalAlignment','Center','VerticalAlignment','Bottom')

plot(primary_1(1),primary_1(2),'b.','Markersize',20)
plot(primary_2(1),primary_2(2),'g.','Markersize',20)

text(primary_1(1),primary_1(2), '$m_1$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',12,'FontName','Times')
text(primary_2(1),primary_2(2), '$m_2$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',12,'FontName','Times')