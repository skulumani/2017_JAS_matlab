%% Plot U3 Section
U3_fig = figure();
set(0,'CurrentFigure',U3_fig)
hold on

for ii = 1:constants.manifold_steps
    uspms = L1_us_manifold_pos_state(:,:,ii);
    usnms = L1_us_manifold_neg_state(:,:,ii);
 
    spms = L1_s_manifold_pos_state(:,:,ii);
    snms = L1_s_manifold_neg_state(:,:,ii);
    
    max_t_usnm = L1_manifold.us_manifold_neg_time_cross(3,ii);
    [~,usnms_cross_index ] = min(abs(L1_us_manifold_neg_time(:,ii)  - max_t_usnm));
    
    plot(usnms(1:usnms_cross_index,1),usnms(1:usnms_cross_index,2),'r')
end

plot_trajectories(t_1, state_1, constants.e_desired, U3_fig, constants)

line([1-constants.mu 1-constants.mu],[0 -0.12],'Linewidth',4,'Color','k')
title('U3 Invariant Manifolds')
% Plot Poincare Section

U3_sect_fig = figure(2);
set(0,'CurrentFigure',U3_sect_fig)
hold on
title('U3 Poincare Section y < 0 x = 1 - \mu')
xlabel('x Axis (nondim)')
ylabel('xdot Axis (nondim)')
grid on
for ii = 1:constants.manifold_steps
    
    
    usnm_cross_state = L1_manifold.us_manifold_neg_state_cross(3,:,ii);
    
    plot(usnm_cross_state(1),usnm_cross_state(3),'r.','Markersize',20)
    
end
