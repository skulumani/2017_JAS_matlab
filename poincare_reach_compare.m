% plot the Poincare section to compare the reachability set against thhose
% of the invariant manifold 
% compare times of flight

% also input the target manifold
function [reach_poincare]=poincare_reach_compare(sol_output,traj_fig, poincare_fig1)

constants = sol_output(1).constants;

% propogate each planar orbit with no control and plot to visualize it
constants.control_switch = 'off';




% parse out the states for each theta angle
num_steps = sol_output(1).constants.num_steps;
num_seg = sol_output(1).constants.num_seg;
num_theta = length(sol_output);
num_states = sol_output(1).constants.num_states;

reach_struct(num_theta) = struct('state',[],'costate',[],'reach_end',[]);
for ii = 1:num_theta % loop over theta angles (poincare directions)
    state = zeros(num_steps,num_states);
    costate = zeros(num_steps,num_states);
    
    % loop over the segments and combine trajectories into a big array
    for jj = 1:num_seg
       x_i = sol_output(ii).x_i;
       h_i = sol_output(ii).h_i;
       start_idx = (jj-1)*num_steps/num_seg+1;
       end_idx = start_idx-1+num_steps/num_seg;
       state(start_idx:end_idx,:) = x_i(:,:,jj);
       costate(start_idx:end_idx,:) = h_i(:,:,jj);
    end
    
    reach_struct(ii).state= state;
    reach_struct(ii).costate = costate;
    reach_struct(ii).reach_end = [state(end,:) costate(end,:)];
    % plot the trajectories on the same plot of the moon periodic orbit
    set(0,'CurrentFigure',traj_fig);
    plot(state(:,1),state(:,2),'.')
    set(0,'CurrentFigure',poincare_fig1);
    plot(state(end,1),state(end,3),'.','Markersize',20)
    
end


reach_poincare = cat(1,reach_struct(:).reach_end);


% %% load the manifold from the periodic orbit
% load './u=05/l1_manifold.mat'
% crossing_sel = 3;
% % plot manifold poincare section and trajectories
% set(0,'CurrentFigure',traj_fig2)
% hold on
% 
% for ii = 1:5:constants.manifold_steps
%     uspms = L1_us_manifold_pos_state(:,:,ii);
%     usnms = L1_us_manifold_neg_state(:,:,ii);
%  
%     spms = L1_s_manifold_pos_state(:,:,ii);
%     snms = L1_s_manifold_neg_state(:,:,ii);
%     
%     max_t_usnm = L1_manifold.us_manifold_neg_time_cross(crossing_sel,ii);
%     [~,usnms_cross_index ] = min(abs(L1_us_manifold_neg_time(:,ii)  - max_t_usnm));
%     
%     plot(usnms(1:usnms_cross_index,1),usnms(1:usnms_cross_index,2),'r')
% end
% 
% plot(state_1(:,1),state_1(:,2),'k','linewidth',4)
% plot_trajectories(t, state, energyconst(moon_x0',constants.mu), traj_fig2, constants)
% line([0.8352 1.1],[0 0],'Linewidth',4,'Color','k')
% % Plot Poincare Section
% 
% 
% set(0,'CurrentFigure',poincare_fig1)
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
% 
% set(0,'CurrentFigure',poincare_fig2)
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
