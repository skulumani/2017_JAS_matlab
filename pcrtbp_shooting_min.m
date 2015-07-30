% 22 July 2015
% multiple shooting to minimize to a desired target
function pcrtbp_shooting_min(min_reach)
% close all
clc
% clear all

% constants and problem setup
constants = min_reach.constants;
constants.optfsolve = optimoptions(@fsolve,'Display','iter-detailed','TolFun',1e-6,'TolX',1e-6,...
    'Algorithm', 'trust-region-reflective','Jacobian','off',...
    'DerivativeCheck','off');
% generate an initial periodic orbit

% iterations from um = 0.1
% min state from l1_reach1.mat (second iteration)
% x0_i=[0.854230866854948   0.000002940915581   0.017850049210949 -0.150017324041792];
% min state from l1_reach_second (third iteration)
% x0_i = [0.929235180302524  -0.001184714718967   0.480459413062010   0.012501857639028];
% periodic orbit

% iterations from um = 0.25
% x0_i = [0.843741045660152   0.000014266827793   0.056397653670216  -0.066280767615112];

x0_i = [0.815614054266804 0 0 0.192227407664904]; % initial condition on periodic orbit
T_i = 1.407478324303006; % half period of initial periodic orbit

% get the h0 from the reach state that is closest to the intersection point
h0 = min_reach.h0;

% beta = [0.1;0.1;0.1;0.1]; % terminal constraint to match a target
num_con = 0;
constants.num_con = num_con;


% generate the sub intervals and a guess of trajectories
num_seg = 4; % should be cleanly divisble by the number of steps
t = min_reach.t;

% load x_i and h_i from reachable set
x_i = min_reach.x_i;
h_i = min_reach.h_i;

h = t(1,2)-t(1,1); % step size
constants.h = h;
constants.t = t;

num_states = size(x_i,2);
constants.num_states = num_states;

num_mid = num_seg-1;

% fixed final state to minimize towards ( from um = 0.5 iteration)
constants.xt = [0.927258145099694 0  -0.036751604638168   -0.369096025781395];
% solver loop
% eps=1e-4;
% max_iter = 50;

constants.x0 = x0_i;

% all the interior patch points for the multiple shooting
xm = squeeze(x_i(end,:,1:end-1));
hm = squeeze(h_i(end,:,1:end-1));

xm = reshape(xm,num_states,num_seg-1);
hm = reshape(hm,num_states,num_seg-1);

%% loop over theta
% solution array

constants.control_switch = 'on';

% minimize to a desired state
% form xg
% [h0, xm, hm, beta]
% xg = zeros(num_states + num_con + num_mid*2*num_states,1);
%     for mid = 1:num_mid
%         srow_idx = num_states+num_states*2*(mid-1)+1;
%         erow_idx8 = srow_idx + 2*num_states - 1;
%         
%         xg(srow_idx:erow_idx8) = [xm(:,mid);hm(:,mid)];
%     end
%     xg(1:num_states) = h0;
%     xg(length(xg)-num_con+1:end) = beta;

% fsolve
% [xg,fval,exitflag,output] = fsolve(@(xg)objective_minimize(xg, x0, constants),xg,constants.optfsolve);
[xg, fval, exitflag, output] = fsolve(@(h0)obj(h0,[0 t(num_seg,end)],constants),h0, constants.optfsolve);

% [x_i, h_i, xm, hm, h0, beta] = prop_seg(xg,x0,constants);
%     plot_seg(t,x_i,h_i,xm,hm, x0, h0);

%     % calculate all the minus states
xminus = squeeze(x_i(end,:,:));
hminus = squeeze(h_i(end,:,:));
% compare to the plus states
diff_x = xminus(:,1:size(xm,2)) - xm;
diff_h = hminus(:,1:size(hm,2)) - hm;
% extract outputs
sol_output.xg_sol = xg;
sol_output.fval = fval;
sol_output.exitflag = exitflag;
sol_output.x_i = x_i;
sol_output.h_i = h_i;
sol_output.xm = xm;
sol_output.hm = hm;
sol_output.h0 = h0;
sol_output.diff_x = diff_x;
sol_output.diff_h = diff_h;
sol_output.costate_necc = hminus(:,end) - ( xminus(:,end) - constants.xt');
sol_output.t = t;
sol_output.constants = constants;
sol_output.x0 = x0;

keyboard
plot_output(sol_output,constants)

% save('l1_reach_025','sol_output')
end

function [con, jac] = objective_minimize(xg, x_0, constants);
% minimize to a desired state
% output the constraint and jacobian
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;
num_steps = constants.num_steps;
num_con = constants.num_con;

t = constants.t;
% xg is a big vector of all the unknowns to be solved for
% unpack xg
eye8x8 = eye(8,8);

h0 = xg(1:num_states);
% beta = xg(length(xg)-(num_con-1):end);
beta = 1;
% parse out xm, hm interior constraint points
xm = zeros(num_states,num_mid);
hm = zeros(num_states,num_mid);

for mid = 1:num_mid
    srow_idx = num_states+num_states*2*(mid-1)+1;
    erow_idx8 = srow_idx + 2*num_states - 1;
    
    mid_comb = xg(srow_idx:erow_idx8);
    xm(:,mid) = mid_comb(1:num_states);
    hm(:,mid) = mid_comb(num_states+1:end);
    
end

% propogate state, costate, sensitivities for each segment
% Generate a new complete trajectory

constants.control_switch = 'on';

% propogate each segment/stage using the linearized system to generate
% x_i, lambda_i
x_i = zeros(num_steps/num_seg,num_states,num_seg); % segment trajectories
h_i = zeros(num_steps/num_seg,num_states,num_seg);

Pxh0 = zeros(num_states);
Phh0 = eye(num_states);
Pxx0 = eye(num_states);
Phx0 = zeros(num_states);

sensf.Pxhf= zeros(num_states,num_states, num_seg); % final segment sensitivities
sensf.Phhf = zeros(num_states,num_states, num_seg);
sensf.Pxxf = zeros(num_states,num_states, num_seg);
sensf.Phxf = zeros(num_states,num_states, num_seg);

sensf.Pxh = zeros(num_states,num_states,num_steps/num_seg,num_seg);
sensf.Phh = zeros(num_states,num_states,num_steps/num_seg,num_seg);
sensf.Pxx = zeros(num_states,num_states,num_steps/num_seg,num_seg);
sensf.Phx = zeros(num_states,num_states,num_steps/num_seg,num_seg);
for seg = 1:num_seg
    % loop over each STM to calculate the state
    if seg == 1
        [~, z,sens,~] = pcrtbp_var_shooting(x_0,h0',Pxh0, Phh0, Pxx0, Phx0, 0,t(1,end), 'trap', num_steps/num_seg, constants);
    else
        [~, z,sens,~] = pcrtbp_var_shooting(xm(:,seg-1),hm(:,seg-1),Pxh0, Phh0, Pxx0, Phx0, t(seg,1),t(seg,end), 'trap', num_steps/num_seg, constants);
    end
    
    % store the propogated state/costate/sensitivites into arrays for each
    % segment
    x_i(:,:,seg) = z(:,1:4);
    h_i(:,:,seg) = z(:,5:8);
    
    sensf.Pxhf(:,:,seg) = sens.Pxh(:,:,end);
    sensf.Phhf(:,:,seg) = sens.Phh(:,:,end);
    sensf.Pxxf(:,:,seg) = sens.Pxx(:,:,end);
    sensf.Phxf(:,:,seg) = sens.Phx(:,:,end);
    
    sensf.Pxh(:,:,:,seg) = sens.Pxh;
    sensf.Phh(:,:,:,seg) = sens.Phh;
    sensf.Pxx(:,:,:,seg) = sens.Pxx;
    sensf.Phx(:,:,:,seg) = sens.Phx;
    
end

% form constraint and jacobian by looping over the stages

con = zeros(num_con+num_states+num_mid*2*num_states,1);
jac = zeros(num_con+num_states+num_mid*2*num_states, num_con+num_states+num_mid*2*num_states);

% loop over the interior segments
for seg = 2:num_seg-1
    % indicies for forming the jacobian and constraint equations
    srow_idx = num_states*2*(seg-1)+1;
    erow_idx = srow_idx + 2*num_states - 1;
    
    scol_idx = num_states+ 2*num_states*(seg-2)+1;
    ecol_idx = scol_idx + 4*num_states - 1;
    
    [con_seg, jac_seg] = seq_newton_minimize(seg, x_i, h_i,sensf, x_0,xm, hm, beta,h0, constants); % compute stage con and jacobian
    % combine into a large g, dgdx vectors for output
    con(srow_idx:erow_idx) = con_seg;
    jac(srow_idx:erow_idx,scol_idx:ecol_idx) = [jac_seg -eye8x8];
end

% compute values for the first stage

[con_seg, jac_seg] = seq_newton_minimize(1, x_i, h_i,sensf, x_0,xm, hm,beta, h0,constants);
con(1:2*num_states) = con_seg;
jac(1:2*num_states,1:3*num_states) = [jac_seg -eye8x8];

% compute for last stage
[con_seg, jac_seg] = seq_newton_minimize(num_seg, x_i, h_i,sensf, x_0,xm, hm,beta, h0,constants);

con(length(con)-(num_states+num_con-1):end) = con_seg;
jac(size(jac,1)-(num_states+num_con-1):end,size(jac,2)-(num_con+2*num_states-1):end) = jac_seg;

end

function [g, dgdx] = seq_newton_minimize(stage, x_i, h_i, sensf, x_0, xm, hm, beta,h0, constants)
% minimize to a target state after reachability sections intersect
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;

% put in logic to determine the number of constraints
num_con = constants.num_con;

% logic if the stage is the first or last stage and to form g, dgdx
if stage == num_seg % last stage
    % calculate the constraint and jacobian of constraint
    xf = x_i(end,:,stage);
    hf = h_i(end,:,stage);
    
    % terminal constraints
    %     [m, dmdx,dm1xdx,dm2xdx] = constraints(xf,hf, constants);
    %     [m, dmdx,dm1xdx,dm2xdx] = constraints_nodenom(xf',hf',constants);
    % newton iteration
    % form the vector g
%     g1 = hf' - beta;
%     g2 =( xf' - constants.xt');
%     
%     
%     g = [g1;g2];
g = xf' - constants.xt';
    %     g = hf + Q *(xf-constants.xcf);
    % form dgdx
%     dg1dxm = sensf.Phxf(:,:,stage) ;
%     dg1dhm = sensf.Phhf(:,:,stage) ;
%     dg1dbeta = -eye(4,4);
%     
%     dg2dxm = sensf.Pxxf(:,:,stage);
%     dg2dhm = sensf.Pxhf(:,:,stage);
%     dg2dbeta = zeros(4,4);
%     
%     dgdx = [dg1dxm dg1dhm dg1dbeta;...
%         dg2dxm dg2dhm dg2dbeta];
    
    dgdx = [ sensf.Pxxf(:,:,stage) sensf.Pxhf(:,:,stage)];
elseif stage == 1 % first stage
    
    % calculate the constraint and jacobian of constraint
    xf = x_i(end,:,stage);
    hf = h_i(end,:,stage);
    
    g1 = xf' - xm(:,stage);
    g2 = hf' - hm(:,stage);
    
    g = [g1;g2];
    dgdx = [sensf.Pxhf(:,:,stage); sensf.Phhf(:,:,stage)];
    
else % all the interior stages
    
    xminus = x_i(end,:,stage);
    hminus = h_i(end,:,stage);
    
    g1 = xminus' - xm(:,stage);
    g2 = hminus' - hm(:,stage);
    
    g = [g1;g2];
    
    dgdx = [sensf.Pxxf(:,:,stage) sensf.Pxhf(:,:,stage);sensf.Phxf(:,:,stage) sensf.Phhf(:,:,stage)];
end

end % end seq newton minimize function

function [x_i, h_i, xm, hm, h0, beta] = prop_seg(xg,x0,constants)
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;
num_steps = constants.num_steps;
num_con = constants.num_con;

t = constants.t;
% xg is a big vector of all the unknowns to be solved for
% unpack xg
eye8x8 = eye(8,8);

h0 = xg(1:num_states);
beta = xg(length(xg)-1:end);
% parse out xm, hm interior constraint points
xm = zeros(num_states,num_mid);
hm = zeros(num_states,num_mid);

for mid = 1:num_mid
    srow_idx = num_states+num_states*2*(mid-1)+1;
    erow_idx8 = srow_idx + 2*num_states - 1;
    
    mid_comb = xg(srow_idx:erow_idx8);
    xm(:,mid) = mid_comb(1:num_states);
    hm(:,mid) = mid_comb(num_states+1:end);
    
end

% propogate state, costate, sensitivities for each segment
% Generate a new complete trajectory

constants.control_switch = 'on';

% propogate each segment/stage using the linearized system to generate
% x_i, lambda_i
x_i = zeros(num_steps/num_seg,4,num_seg);
h_i = zeros(num_steps/num_seg,4,num_seg);
for seg = 1:num_seg
    % loop over each STM to calculate the state
    state = zeros(size(x_i(:,:,seg)));
    costate = zeros(size(h_i(:,:,seg)));
    if seg == 1
        
        [~, z,~] = pcrtbp_optimal_var(x0,h0',0,t(1,end), 'trap', num_steps/num_seg, constants);
    else
        
        [~, z,~] = pcrtbp_optimal_var(xm(:,seg-1),hm(:,seg-1),t(seg,1),t(seg,end), 'trap', num_steps/num_seg, constants);
    end
    
    
    x_i(:,:,seg) = z(:,1:4);
    h_i(:,:,seg) = z(:,5:8);
    
    
end

end

function plot_output(sol_output,constants)
num_steps = sol_output.constants.num_steps;
num_seg = sol_output.constants.num_seg;

num_states = sol_output.constants.num_states;

figure(1)
grid on
hold all
    state = zeros(num_steps,num_states);
    costate = zeros(num_steps,num_states);
% loop over the segments and combine trajectories into a big array
    for ii = 1:constants.num_seg
       x_i = sol_output.x_i;
       h_i = sol_output.h_i;
       start_idx = (ii-1)*num_steps/num_seg+1;
       end_idx = start_idx-1+num_steps/num_seg;
       state(start_idx:end_idx,:) = x_i(:,:,ii);
       costate(start_idx:end_idx,:) = h_i(:,:,ii);
    end
    plot(state(:,1),state(:,2),'g')
end

function F = obj(h0,tspan,constants)


% integrate the equations of motion for the given h_0 
constants.control_switch = 'on';
[~, state,~] = pcrtbp_optimal_var(constants.x0,h0',0,tspan(end), 'trap', 1000, constants);
% subtract the final state from the desired final state 
F = ( (state(end,1:4)'- constants.xt') );


end