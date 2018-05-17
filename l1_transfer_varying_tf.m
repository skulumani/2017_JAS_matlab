
%% vary the maximum control to see effect on reachability set

clear all
close all
clc

um_array = [0.5, 0.6, 0.7, 0.8, 0.9];
color = ['r', 'g', 'b', 'y', 'm'];
marker = ['o', '+', '*', 'x', 's'];

initial_condition = [  0.815614054266804, 0, 0, 0.192227407664904];
reach_time = 1.307478324303006;
reach_time_array = 1.307478324303006 * [0.95, 1.0, 1.05, 1.1];

%% now loop over different terminal times

for ii = 1:size(reach_time_array, 2)
    sol_output = pcrtbp_shooting(initial_condition, reach_time_array(ii), 0.25);
    save(['./data/l1_varying_tf_um_25/l1_reach_', num2str(reach_time_array(ii)), '.mat'], 'sol_output');
end

fprintf('Done with tf loop\n')
