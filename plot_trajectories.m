% function to plot a state trajectory

% created 12 september
% 15 september added energy level input and figure handle

% 24 october 2014 - modifications to plot trajectories intelligently based
% on number of states (pcrtbp vs crtbp)
function plot_trajectories(t, state, E, fig_handle, constants)

switch constants.plot_scale
    case 'dim'
        l_scale = constants.l_scale;
        units_label = '(km)';
    case 'nondim'
        l_scale = 1;
        units_label = '(nondim)';
end

mu = constants.mu;

grid_spacing = 0.01;

L_points = libration_points(mu);

L1 = [L_points(1,:) 0];
L2 = [L_points(2,:) 0];
L3 = [L_points(3,:) 0];
L4 = [L_points(4,:) 0];
L5 = [L_points(5,:) 0];

% % energies of each libartion point
% [C1 E1] = energyconst([L1 zeros(1,3)],mu);
% [C2 E2] = energyconst([L2 zeros(1,3)],mu);
% [C3 E3] = energyconst([L3 zeros(1,3)],mu);
% [C4 E4] = energyconst([L4 zeros(1,3)],mu);
% [C5 E5] = energyconst([L5 zeros(1,3)],mu);


% E = E2 + 0.1*abs(E2-E3);
% E=constants.e_desired;

%Contour plot of the zero velocity energy curve
[X,Y]=meshgrid(-2:grid_spacing:2);
r1=((X + mu).^2+Y.^2 ).^(1/2);
r2=((X + mu-1).^2+Y.^2 ).^(1/2);
C_grid=2*((1/2)*(X.^2 + Y.^2)+(1-mu)./r1+mu./r2); % jacobi constant
E_grid = C_grid/-2;

primary_1 = l_scale*[-mu 0 0];
primary_2 = l_scale* [1-mu 0 0];

set(0,'CurrentFigure',fig_handle)
hold all
grid on
title_string = ['Trajectory ($\mu = $' num2str(mu) ')'];
title(title_string,'interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
xlabel(['X' units_label],'interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
ylabel(['Y' units_label],'interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
zlabel(['Z' units_label],'interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')

% contour(X,Y,E_grid,[E,E],'r')

plot3(l_scale*L1(1), 0,0, 'k.', l_scale*L2(1), 0,0, 'k.', l_scale*L3(1), 0,0,'k.','Markersize',20)
plot3(l_scale*L4(1), l_scale*L4(2),0,'k.','Markersize',20)
plot3(l_scale*L5(1), l_scale*L5(2),0,'k.','Markersize',20)

text(l_scale*L1(1),0, '$L_1$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
text(l_scale*L2(1),0, '$L_2$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
text(l_scale*L3(1),0, '$L_3$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
text(l_scale*L4(1),l_scale*L4(2), '$L_4$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
text(l_scale*L5(1),l_scale*L5(2), '$L_5$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')

plot(primary_1(1),primary_1(2),'b.','Markersize',20)
plot(primary_2(1),primary_2(2),'b.','Markersize',20)

text(primary_1(1),primary_1(2), '$m_1$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')
text(primary_2(1),primary_2(2), '$m_2$','HorizontalAlignment','Center','VerticalAlignment','Bottom','interpreter','latex','FontUnits','points','FontSize',9,'FontName','Times')

% switch to plot trajectory
num_states = size(state,2);

switch num_states
    case 4
        plot(l_scale*state(:,1),l_scale*state(:,2))
    case 6
        plot3(l_scale*state(:,1),l_scale*state(:,2),l_scale*state(:,3))
end

% plot(l_scale*-mu,0,'r*',l_scale*(1-mu),0,'r*')
