% 26 October 2015
% parse out a stable manifold

function manifold_poincare = manifold_parse(traj_fig, poincare_fig)
load './l1_manifold_geo.mat'

% number of crossings of the x axis to plot
crossing_sel = 3;

% store all manifold poincare section states to compare with the
% reachability set
manifold_poincare = zeros(constants.manifold_steps,4);
% plot the stable manifold of L1
for ii = 1:1:constants.manifold_steps
    uspms = L1_us_manifold_pos_state(:,:,ii);
    usnms = L1_us_manifold_neg_state(:,:,ii);
    
    spms = L1_s_manifold_pos_state(:,:,ii);
    snms = L1_s_manifold_neg_state(:,:,ii);
    
    % maximum time along the manifold (number of x axis crossings)
    %     max_t_uspm = L1_manifold.us_manifold_pos_time_cross(crossing_sel,ii);
    %     [~,uspms_cross_index ] = min(abs(L1_us_manifold_pos_time(:,ii)  - max_t_uspm));
    %
    %     max_t_usnm = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
    %     [~,usnms_cross_index ] = min(abs(L1_us_manifold_neg_time(:,ii)  - max_t_usnm));
    %
    %     max_t_spm = L1_manifold.s_manifold_pos_time_cross(crossing_sel,ii);
    %     [~,spms_cross_index ] = min(abs(L1_s_manifold_pos_time(:,ii)  - max_t_spm));
    %
    %     max_t_snm = L1_manifold.s_manifold_neg_time_cross(crossing_sel,ii);
    %     [~,snms_cross_index ] = min(abs(L1_s_manifold_neg_time(:,ii)  - max_t_snm));
    
    % find index of crossing state that hits the U1 section y = 0 x < 0
    snm_index = find(L1_manifold.s_manifold_neg_state_cross(:,1,ii) < 0);
    max_t_snm = L1_manifold.s_manifold_neg_time_cross(snm_index(1),ii);
    [~,snms_cross_index ] = min(abs(L1_s_manifold_neg_time(:,ii)  - max_t_snm));
    
    snm_cross_state = L1_manifold.s_manifold_neg_state_cross(snm_index(1),:,ii);
    manifold_poincare(ii,:) = snm_cross_state;
    % plot the manifold
    %     plot(uspms(1:uspms_cross_index,1),uspms(1:uspms_cross_index,2),'r')
    %     plot(usnms(1:usnms_cross_index,1),usnms(1:usnms_cross_index,2),'r')
    
    %     plot(spms(1:spms_cross_index,1),spms(1:spms_cross_index,2),'g')
    set(0,'CurrentFigure',traj_fig)
    plot(snms(1:snms_cross_index,1),snms(1:snms_cross_index,2),'g')
    set(0,'CurrentFigure',poincare_fig)
    plot(snm_cross_state(1),snm_cross_state(3),'g.','Markersize',20)
end

% plot(state_1(:,1),state_1(:,2),'k','linewidth',4)

% line([0.8352 1.1],[0 0],'Linewidth',4,'Color','k')
% Plot Poincare Section


%
%
% for ii = 1:5:constants.manifold_steps
%
%
%     usnm_cross_time = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
%     usnm_cross_state = L1_manifold.us_manifold_neg_state_cross(crossing_sel,:,ii);
%
%     plot(usnm_cross_state(1),usnm_cross_state(3),'g.','Markersize',20)
%     text(usnm_cross_state(1),usnm_cross_state(3), num2str(usnm_cross_time))
% end