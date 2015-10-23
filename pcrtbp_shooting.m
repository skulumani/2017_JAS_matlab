% 8 July 2015
% Multiple shooting for PCRTBP Variational problem. Really hope this works
% finally
function [sol_output]= pcrtbp_shooting(initial_condition, reach_time)
% close all
% clc
% clear all

% constants and problem setup
constants = crtbp_constants;
constants.optfsolve = optimoptions(@fsolve,'Display','iter','TolFun',1e-4,'TolX',1e-4,...
    'MaxIter',5000,'MaxFunEvals',5000, 'Algorithm', 'trust-region-dogleg','Jacobian','off',...
    'DerivativeCheck','off');
constants.um = 0.5;

num_steps = 1000;
constants.num_steps = num_steps;

% [L_points, ~] = libration_points(constants.mu);
% constants.center_vec = L_points(1,:);
constants.center_vec = [-constants.mu 0]; % centered at the Earth

constants.alpha_d = 0*pi/180 ; % angle that defines the poincare section

x0_i = initial_condition;
T_i = reach_time;

% propogate the intial and target orbits to the poincare section using
% variational integrator
constants.control_switch = 'off';
h0_i = 0.1*ones(1,4);

[t_i, state_i,~] = pcrtbp_optimal_var(x0_i,h0_i,0,T_i, 'trap', num_steps, constants);


% generate initial guess of state and costate histories
% from periodic orbit
h01 = [   0.032389612643357;...
   0.010508161299126;...
   0.000547425479709;...
   0.001064218889245]; 

beta1 = [   1.256812564325579;...
  -0.503978575859051];

num_con = 2;
constants.num_con = num_con;

% initial guess of state and costate histories
constants.control_switch = 'off';
tf = T_i;

% generate the sub intervals and a guess of trajectories
num_seg = 4; % should be cleanly divisble by the number of steps
t = zeros(num_seg, num_steps/num_seg);
x_i = zeros(num_steps/num_seg,4,num_seg);
h_i = zeros(num_steps/num_seg,4,num_seg);

constants.num_seg = num_seg;

for ii = 1:num_seg
    t(ii,:) = linspace((ii-1)*tf/num_seg,ii*tf/num_seg,num_steps/num_seg);
    % propogate the no control trajectory to generate the interior guess points
    if ii == 1
        % propogate from x0, h0
        [~, state,~] = pcrtbp_optimal_var(x0_i,h0_i,0,t(ii,end), 'trap', num_steps/num_seg, constants);
        x_i(:,:,ii) = state(:,1:4);
        h_i(:,:,ii) = state(:,5:8);
    else
        [~, state,~] = pcrtbp_optimal_var(x_i(end,:,ii-1),h_i(end,:,ii-1),t(ii-1,end),t(ii,end), 'trap', num_steps/num_seg, constants);
        x_i(:,:,ii) = state(:,1:4);
        h_i(:,:,ii) = state(:,5:8);
    end
end
h = t(1,2)-t(1,1); % step size
constants.h = h;
constants.t = t;

num_states = size(x_i,2);
constants.num_states = num_states;

num_mid = num_seg-1;
Q = [1 0 0 0;0 0 0 0;0 0 1 0;0 0 0 0];
% no control final state to maximize against
constants.xcf = state_i(end,1:4)';

% solver loop
% eps=1e-4;
% max_iter = 50;

x0 = x0_i;

% all the interior patch points for the multiple shooting
xm = squeeze(x_i(end,:,1:end-1));
hm = squeeze(h_i(end,:,1:end-1));

xm = reshape(xm,num_states,num_seg-1);
hm = reshape(hm,num_states,num_seg-1);

%% loop over theta
% solution array

theta =linspace(0,360,180)*pi/180;

constants.control_switch = 'on';
sol_output(length(theta)) = struct('theta',[],'xg_sol',[], 'fval',[],'exitflag',[],...
    'x_i',[],'h_i',[],'xm',[],'hm',[],'h0',[],'beta',[],'diff_x',[],...
    'diff_h',[],'m',[],'costate_necc',[],'t',[],'constants',[],'x0',[]);


for theta_ind = 1:length(theta)% loop over theta angles
 constants.theta_d = theta(theta_ind);
    % loop or fsolve
    
    % form the xg guess vector
    % [h0, xm, hm, beta]
    xg = zeros(num_states + num_con + num_mid*2*num_states,1);
    for mid = 1:num_mid
        srow_idx = num_states+num_states*2*(mid-1)+1;
        erow_idx8 = srow_idx + 2*num_states - 1;
        
        xg(srow_idx:erow_idx8) = [xm(:,mid);hm(:,mid)];
    end
    xg(1:num_states) = h01;
    xg(length(xg)-num_con+1:end) = beta1;
    % iterate to solve for xg
    % test the jacobian
%     jacobian_test(xg,Q, x0,constants)
%     % newton iteration loop here
%     iter = 1;
%     g = 1;
%     while norm(g) > 1e-3 && iter < 100
%         % calculate g, dgdx and use a newton iteration to solve
%         
%         [g, dgdx] = objective(xg, Q, x0, constants);
%         
%         del = -0.5 * (dgdx \ g);
%         xg = xg + del;
%         clc
%         fprintf('del norm %5.4f\n',norm(del));
%         fprintf('g norm %5.4f\n',norm(g));
%         iter = iter + 1;
%     end

    [xg,fval,exitflag,output] = fsolve(@(xg)objective(xg,Q, x0, constants),xg,constants.optfsolve);
    
    [x_i, h_i, xm, hm, h0, beta] = prop_seg(xg,x0,constants);
    %     plot_seg(t,x_i,h_i,xm,hm, x0, h0);
    
    % calculate all the minus states
    xminus = squeeze(x_i(end,:,:));
    hminus = squeeze(h_i(end,:,:));
    % compare to the plus states
    diff_x = xminus(:,1:size(xm,2)) - xm;
    diff_h = hminus(:,1:size(hm,2)) - hm;
    
    % neccessary condition on terminal costate
    [m, dmdx,~,~] = constraints_nodenom(xminus(:,end),hminus(:,end),constants);
    
    costate_necc = hminus(:,end) + Q * ( xminus(:,end) - constants.xcf) - dmdx' * beta;
    
    % save to output structure for later plotting
    sol_output(theta_ind).theta = constants.theta_d;
    sol_output(theta_ind).xg_sol = xg;
    sol_output(theta_ind).fval = fval;
    sol_output(theta_ind).exitflag = exitflag;
    sol_output(theta_ind).x_i = x_i;
    sol_output(theta_ind).h_i = h_i;
    sol_output(theta_ind).xm = xm;
    sol_output(theta_ind).hm = hm;
    sol_output(theta_ind).h0 = h0;
    sol_output(theta_ind).beta = beta;
    sol_output(theta_ind).diff_x = diff_x;
    sol_output(theta_ind).diff_h = diff_h;
    sol_output(theta_ind).m = m;
    sol_output(theta_ind).costate_necc = costate_necc;
    sol_output(theta_ind).t = t;
    sol_output(theta_ind).x0 = x0;
    fprintf('Theta %5.4f \n',constants.theta_d);
end % theta angle loop

keyboard
plot_output(sol_output,constants)

sol_output(1).constants = constants;
% save('l1_reach_025','sol_output')
end

function [con, jac] = objective(xg, Q, x_0,constants)
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
    
    [con_seg, jac_seg] = seq_newton(seg, x_i, h_i,sensf,Q, x_0,xm, hm,beta, h0, constants); % compute stage con and jacobian
    % combine into a large g, dgdx vectors for output
    con(srow_idx:erow_idx) = con_seg;
    jac(srow_idx:erow_idx,scol_idx:ecol_idx) = [jac_seg -eye8x8];
end

% compute values for the first stage

[con_seg, jac_seg] = seq_newton(1, x_i, h_i,sensf,Q, x_0,xm, hm, beta, h0,constants);
con(1:2*num_states) = con_seg;
jac(1:2*num_states,1:3*num_states) = [jac_seg -eye8x8];

% compute for last stage
[con_seg, jac_seg] = seq_newton(num_seg, x_i, h_i,sensf,Q, x_0,xm, hm, beta, h0,constants);

con(length(con)-(num_states+num_con-1):end) = con_seg;
jac(size(jac,1)-(num_states+num_con-1):end,size(jac,2)-(num_con+2*num_states-1):end) = jac_seg;


end % end of objective function




function [g, dgdx] = seq_newton(stage, x_i, h_i, sensf,Q, x_0,xm, hm, beta, h0,constants)
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
    [m, dmdx,dm1xdx,dm2xdx] = constraints_nodenom(xf',hf',constants);
    % newton iteration
    % form the vector g
    g1 =  hf' + Q * ( xf' - constants.xcf) - dmdx' * beta;
    g2 = m;
    
    g = [g1;g2];
    %     g = hf + Q *(xf-constants.xcf);
    % form dgdx
    dg1dxm = sensf.Phxf(:,:,stage) + Q*sensf.Pxxf(:,:,stage) - dm1xdx*sensf.Pxxf(:,:,stage)*beta(1) - dm2xdx*sensf.Pxxf(:,:,stage)*beta(2);
    dg1dhm = sensf.Phhf(:,:,stage) + Q*sensf.Pxhf(:,:,stage) - dm1xdx*sensf.Pxhf(:,:,stage)*beta(1) - dm2xdx*sensf.Pxhf(:,:,stage)*beta(2);
    dg1dbeta = -dmdx';
    
    dg2dxm = dmdx*sensf.Pxxf(:,:,stage);
    dg2dhm = dmdx*sensf.Pxhf(:,:,stage);
    dg2dbeta = zeros(2,2);
    
    dgdx = [dg1dxm dg1dhm dg1dbeta;...
        dg2dxm dg2dhm dg2dbeta];
    
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

end % seq newton function

function jacobian_test(xg,Q, x0,constants)
% test the boundary jacobians
deltax = zeros(size(xg));
deltax(5:8) = 1e-3*rand(4,1);

[g, dgdx] = objective(xg, Q, x0, constants);

% perturb the xm and hm


[g_del, ~] = objective(xg+deltax,  Q, x0, constants);

deltag = g_del-g;
jac = dgdx*deltax;

diffg = abs(deltag-jac);

disp(diffg)
end

function plot_seg(t,x_i,h_i,xm,hm, x0, h0)
% plot the segments to see if they're correct
traj_fig = figure(1);
grid on
hold all

% costate_fig = figure(2);
% grid on
% hold all
% state_fig = figure(3);
% grid on
% hold all
for ii = 1:size(x_i,3)
    set(0,'CurrentFigure',traj_fig);
    plot(x_i(:,1,ii),x_i(:,2,ii))
    if ii == 1
        plot(x0(1),x0(2),'ko')
        text(x0(1),x0(2),sprintf('%s%d%s','x',ii-1,'+'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        statef = [x_i(end,:,ii) h_i(end,:,ii)];
        
    else
        plot(xm(1,ii-1),xm(2,ii-1),'ko')
        text(xm(1,ii-1),xm(2,ii-1),sprintf('%s%d%s','x',ii-1,'+'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        statef = [x_i(end,:,ii) h_i(end,:,ii)];
        
    end
    plot(statef(1),statef(2),'rs','markersize',10)
    text(statef(1),statef(2),sprintf('%s%d%s','x',ii,'-'),'HorizontalAlignment','right','verticalalignment','bottom','interpreter','latex')
    % calculate the x_minus, h_minus terms
    
%     % costate
%     set(0,'CurrentFigure',costate_fig);
%     for ind = 1:size(h_i,2)
%         subplot(4,1,ind)
%         grid on
%         hold all
%         plot(t(ii,:),h_i(:,ind,ii))
%         title(sprintf('%s %d','lambda',ind));
%         
%         if ii == 1
%             plot(t(ii,1),h0(ind),'ko')
%             text(t(ii,1),h0(ind),sprintf('%s%d%s%d%s','$\lambda_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
%         else
%             plot(t(ii,1),hm(ind,ii-1),'ko')
%             text(t(ii,1),hm(ind,ii-1),sprintf('%s%d%s%d%s','$\lambda_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
%         end
%         
%         plot(t(ii,end),statef(4+ind),'rs','markersize',10)
%         text(t(ii,end),statef(4+ind),sprintf('%s%d%s%d%s','$\lambda_{',ind,'}^{',ii,'}-$'),'HorizontalAlignment','right','verticalalignment','bottom','interpreter','latex')
%     end
%     % state
%     set(0,'CurrentFigure',state_fig);
%     for ind = 1:size(h_i,2)
%         subplot(4,1,ind)
%         grid on
%         hold all
%         plot(t(ii,:),x_i(:,ind,ii))
%         title(sprintf('%s %d','x',ind));
%         
%         
%         if ii == 1
%             plot(t(ii,1),x0(ind),'ko')
%             text(t(ii,1),x0(ind),sprintf('%s%d%s%d%s','$x_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
%         else
%             plot(t(ii,1),xm(ind,ii-1),'ko')
%             text(t(ii,1),xm(ind,ii-1),sprintf('%s%d%s%d%s','$x_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
%         end
%         plot(t(ii,end),statef(ind),'rs','markersize',10)
%         text(t(ii,end),statef(ind),sprintf('%s%d%s%d%s','$x_{',ind,'}^{',ii,'}-$'),'HorizontalAlignment','right','verticalalignment','bottom','interpreter','latex')
%     end
    
    
end

% keyboard
end % end of plot segment function

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
poincare_fig = figure(2);
grid on 
hold all

    for ii = 1:length(sol_output)
%         plot_seg(sol_output(ii).t,sol_output(ii).x_i,sol_output(ii).h_i,sol_output(ii).xm,sol_output(ii).hm, sol_output(ii).x0, sol_output(ii).h0);
        set(0,'CurrentFigure',poincare_fig);
        plot(sol_output(ii).x_i(end,1,end),sol_output(ii).x_i(end,3,end),'ro')
    end
end