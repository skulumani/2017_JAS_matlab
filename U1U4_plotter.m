%% Plot U1 and U4 Sections
U1_fig = figure();
set(0,'CurrentFigure',U1_fig)
hold on

for ii = 1:constants.manifold_steps
    uspms = L1_us_manifold_pos_state(:,:,ii);
    usnms = L1_us_manifold_neg_state(:,:,ii);
 
    spms = L1_s_manifold_pos_state(:,:,ii);
    snms = L1_s_manifold_neg_state(:,:,ii);
    
    % find index of crossing state that hits the U1 section y = 0 x < 0
    snm_index = find(L1_manifold.s_manifold_neg_state_cross(:,1,ii) < 0);
    max_t_snm = L1_manifold.s_manifold_neg_time_cross(snm_index(1),ii);
    [~,snms_cross_index ] = min(abs(L1_s_manifold_neg_time(:,ii)  - max_t_snm));
    
    uspm_index = find(L1_manifold.us_manifold_pos_state_cross(:,1,ii) < 0);
    max_t_uspm = L1_manifold.us_manifold_pos_time_cross(uspm_index(1),ii);
    [~,uspms_cross_index ] = min(abs(L1_us_manifold_pos_time(:,ii)  - max_t_uspm));
    
    plot(snms(1:snms_cross_index,1),snms(1:snms_cross_index,2),'g')
    plot(uspms(1:uspms_cross_index,1),uspms(1:uspms_cross_index,2),'r')
end

plot_trajectories(t_1, state_1, constants.e_desired, U1_fig, constants)
plot(state_2(:,1),state_2(:,2),'k')
line([-constants.mu -0.8],[0 0],'Linewidth',4,'Color','k')
title('U1 Invariant Manifolds')
% Plot the U1 Poincare Section

U1_sect_fig = figure();
set(0,'CurrentFigure',U1_sect_fig)
hold on
title('U1 Poincare Section x < 0 y = 0')
xlabel('x Axis (nondim)')
ylabel('xdot Axis (nondim)')
grid on
for ii = 1:constants.manifold_steps
    
    snm_index = find(L1_manifold.s_manifold_neg_state_cross(:,1,ii) < 0);
    uspm_index = find(L1_manifold.us_manifold_pos_state_cross(:,1,ii) < 0);
    snm_cross_state = L1_manifold.s_manifold_neg_state_cross(snm_index(1),:,ii);
    uspm_cross_state = L1_manifold.us_manifold_pos_state_cross(uspm_index(1),:,ii);
    
    plot(uspm_cross_state(1),uspm_cross_state(3),'r.','Markersize',20)
    plot(snm_cross_state(1),snm_cross_state(3),'g.','Markersize',20)
end

% Now plot the U4 exterior section from L2
U4_fig = figure();
set(0,'CurrentFigure',U4_fig)
hold on

for ii = 1:constants.manifold_steps
    uspms = L2_us_manifold_pos_state(:,:,ii);
    usnms = L2_us_manifold_neg_state(:,:,ii);
 
    spms = L2_s_manifold_pos_state(:,:,ii);
    snms = L2_s_manifold_neg_state(:,:,ii);
    
    % find index of crossing state that hits the U4 section y = 0 x < -1
    spm_index = find(L2_manifold.s_manifold_pos_state_cross(:,1,ii) < -1);
    max_t_spm = L2_manifold.s_manifold_pos_time_cross(spm_index(1),ii);
    [~,spms_cross_index ] = min(abs(L2_s_manifold_pos_time(:,ii)  - max_t_spm));
    
    uspm_index = find(L2_manifold.us_manifold_pos_state_cross(:,1,ii) < -1);
    max_t_uspm = L2_manifold.us_manifold_pos_time_cross(uspm_index(1),ii);
    [~,uspms_cross_index ] = min(abs(L2_us_manifold_pos_time(:,ii)  - max_t_uspm));
    
    plot(spms(1:spms_cross_index,1),spms(1:spms_cross_index,2),'g')
    plot(uspms(1:uspms_cross_index,1),uspms(1:uspms_cross_index,2),'r')
end

plot_trajectories(t_1, state_1, constants.e_desired, U4_fig, constants)
plot(state_2(:,1),state_2(:,2),'k')
line([-1 -3.6],[0 0],'Linewidth',4,'Color','k')
title('U4 Invariant Manifolds')
% Plot the U4 Poincare Section

U4_sect_fig = figure();
set(0,'CurrentFigure',U4_sect_fig)
hold on
title('U4 Poincare Section x < -1 y = 0')
xlabel('x Axis (nondim)')
ylabel('xdot Axis (nondim)')
grid on
for ii = 1:constants.manifold_steps
    
    spm_index = find(L2_manifold.s_manifold_pos_state_cross(:,1,ii) < -1);
    uspm_index = find(L2_manifold.us_manifold_pos_state_cross(:,1,ii) < -1);
    spm_cross_state = L2_manifold.s_manifold_pos_state_cross(spm_index(1),:,ii);
    uspm_cross_state = L2_manifold.us_manifold_pos_state_cross(uspm_index(1),:,ii);
    
    plot(uspm_cross_state(1),uspm_cross_state(3),'r.','Markersize',20)
    plot(spm_cross_state(1),spm_cross_state(3),'g.','Markersize',20)
end