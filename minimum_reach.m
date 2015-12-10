% 26 october 2015
% find the minimum of a reachability set to a target set

function [min_reach, min_man, min_traj,min_costate,min_control, min_time,reach_struct] = minimum_reach(sol_output,manifold_poincare, reach_switch)


constants = sol_output(1).constants;

% propogate each planar orbit with no control and plot to visualize it
% constants.control_switch = 'off';

% parse out the states for each theta angle
num_steps = sol_output(1).constants.num_steps;
num_seg = sol_output(1).constants.num_seg;
num_theta = length(sol_output);
num_states = sol_output(1).constants.num_states;
um = constants.um;

reach_struct(num_theta) = struct('state',[],'costate',[],'reach_end',[], 'control',[],'time',[]);
for ii = 1:num_theta % loop over theta angles (poincare directions)
    state = zeros(num_steps,num_states);
    costate = zeros(num_steps,num_states);
    time = zeros(num_steps,1);
    % loop over the segments and combine trajectories into a big array
    for jj = 1:num_seg
        x_i = sol_output(ii).x_i;
        h_i = sol_output(ii).h_i;
        t_i = sol_output(ii).t;
        start_idx = (jj-1)*num_steps/num_seg+1;
        end_idx = start_idx-1+num_steps/num_seg;
        state(start_idx:end_idx,:) = x_i(:,:,jj);
        costate(start_idx:end_idx,:) = h_i(:,:,jj);
        time(start_idx:end_idx,1) = t_i(jj,:)';
    end
    
    reach_struct(ii).state= state;
    reach_struct(ii).costate = costate;
    reach_struct(ii).reach_end = [state(end,:) costate(end,:)];
    
    reach_struct(ii).control = um*costate(:,3:4)./repmat(sqrt(costate(:,3).^2+costate(:,4).^2),1,2);
    reach_struct(ii).time = time;
end


reach_poincare = cat(1,reach_struct(:).reach_end);


% determine the state closest to the target periodic orbit

switch reach_switch
    case 'min_dist'
        % loop over each reach state
        dist_mat = zeros(length(reach_poincare),length(manifold_poincare));
        for ii = 1:length(reach_poincare)
            for jj = 1:length(manifold_poincare)
                % compute the distance (x vs xdot) to the current manifold state
                reach_x = reach_poincare(ii,1);
                reach_xd = reach_poincare(ii,3);
                
                man_x = manifold_poincare(jj,1);
                man_xd = manifold_poincare(jj,3);
                
                dist = sqrt((reach_x-man_x)^2 + (reach_xd-man_xd)^2);
                dist_mat(ii,jj) = dist;
            end
        end
        
        % find the minimum and the row/column of the minimum distance
        [min_row, row_ind] = min(dist_mat);
        [min_col, col_ind] = min(min_row);
        
        % pull out the minimum reach state and corresponding minimum manifold
        min_reach = reach_poincare(row_ind(col_ind),:);
        min_man = manifold_poincare(col_ind,:);
        min_traj = [reach_struct(row_ind(col_ind)).state ];
        min_costate = reach_struct(row_ind(col_ind)).costate;
        min_control = reach_struct(row_ind(col_ind)).control;
        min_time = reach_struct(row_ind(col_ind)).time;
    case 'min_x'
        
        % find the minimum (furthest to the left) final x position of the trajectory
        
        [min_x,x_ind] = min(reach_poincare(:,1));
        dist = zeros(length(manifold_poincare),1);
        % compute the distance from this state to the target manifold
        for ii = 1:length(manifold_poincare)
            reach_x = reach_poincare(x_ind,1);
            reach_xd = reach_poincare(x_ind,3);
            
            man_x = manifold_poincare(ii,1);
            man_xd = manifold_poincare(ii,3);
            
            dist(ii) = sqrt((reach_x-man_x)^2 + (reach_xd-man_xd)^2);
            
        end
        
        [min_row, row_ind] = min(dist);
        min_reach = reach_poincare(x_ind,:);
        min_man = manifold_poincare(row_ind,:);
        min_traj = [reach_struct(x_ind).state ];
        min_costate = reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
    case 'max_x'
        
        % find the maximum (furthest to the right) final x position of the trajectory
        
        [max_x,x_ind] = max(reach_poincare(:,1));
        dist = zeros(length(manifold_poincare),1);
        % compute the distance from this state to the target manifold
        for ii = 1:length(manifold_poincare)
            reach_x = reach_poincare(x_ind,1);
            reach_xd = reach_poincare(x_ind,3);
            
            man_x = manifold_poincare(ii,1);
            man_xd = manifold_poincare(ii,3);
            
            dist(ii) = sqrt((reach_x-man_x)^2 + (reach_xd-man_xd)^2);
            
        end
        
        [min_row, row_ind] = min(dist);
        min_reach = reach_poincare(x_ind,:);
        min_man = manifold_poincare(row_ind,:);
        min_traj = [reach_struct(x_ind).state];
        min_costate =  reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
        
end

