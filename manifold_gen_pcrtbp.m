function [us_manifold_pos_state, us_manifold_neg_state, s_manifold_pos_state, s_manifold_neg_state, ...
    us_manifold_pos_time, us_manifold_neg_time, s_manifold_pos_time, s_manifold_neg_time, manifold]=...
    manifold_gen_pcrtbp(per_orbit_x0,per_orbit_T,per_orbit_E,constants)

% september 4 2014
% september 8 2014 error in backwards propogation for stable manifold

% extract some parameters from the constants structure
manifold_steps = constants.manifold_steps;
numsteps = constants.numsteps;
epsilon = constants.epsilon;
manifold_tf_multipler = constants.manifold_tf_multipler;

manifold_tf = manifold_tf_multipler*per_orbit_T;

mu = constants.mu;

ode_model = 'pcrtbp_stm';
[~, state] = trajectory_simulate(per_orbit_x0, [0 per_orbit_T],ode_model,constants);

cross_fcn = @(t, state)manifold_event(t, state,constants);
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',cross_fcn);


% generate the invariant manifolds
% monodromy matrix based on the STM at the half period or simply integrate
% the whole STM for the perioid
switch ode_model
    case 'pcrtbp_stm'
        stm = state(:,5:end);
        state = state(:,1:4);
        
        M = reshape(stm(end,:),4,4);
    otherwise
        A = diag([1 -1 -1 1]);
        M = A/phi_cross_out*A*phi_cross_out;
end

% eigenvalues of the monodromy matrix
[V, D] = eig(M);

% get stable/unstable eigenvectors at time T
eig_val = diag(D);
[~, us_index] = max(eig_val);
[~, s_index] = min(eig_val);

us_eig_vec0 = V(:,us_index);
s_eig_vec0 = V(:,s_index);

% number of points along periodic orbit to use for manifold
mani_density=round((numsteps-1)/manifold_steps);

us_manifold_pos_initial_state = zeros(manifold_steps,size(state,2));
s_manifold_pos_initial_state = zeros(manifold_steps,size(state,2));

us_manifold_neg_initial_state = zeros(manifold_steps,size(state,2));
s_manifold_neg_initial_state = zeros(manifold_steps,size(state,2));

periodic_orbit_initial = zeros(manifold_steps,size(state,2));
% get the manifold initial conditions
for ii = 1:manifold_steps
    jj = (ii-1)*mani_density+1;
    
    % state on periodic orbit
    initial_xt = state(jj,:);
    % stm at current point on periodic orbit
    stm_xt = reshape(stm(jj,:),4,4);
    % use the STM to move the initial eigenvector all over the periodic orbit
    us_eig_vec = stm_xt*us_eig_vec0; 
    s_eig_vec = stm_xt*s_eig_vec0;
    
    % renormalize
    us_eig_vec = us_eig_vec/norm(us_eig_vec);
    s_eig_vec = s_eig_vec/norm(s_eig_vec);
    % generate the initial perturbations for the manifolds
    us_manifold_pos_initial_state(ii,:) =  initial_xt + epsilon*us_eig_vec';
    us_manifold_neg_initial_state(ii,:) =  initial_xt - epsilon*us_eig_vec';
    
    s_manifold_pos_initial_state(ii,:) =  initial_xt + epsilon*s_eig_vec';
    s_manifold_neg_initial_state(ii,:) =  initial_xt - epsilon*s_eig_vec';
    
    periodic_orbit_initial(ii,:) = initial_xt;
end

us_manifold_pos_state = zeros(numsteps,size(state,2),manifold_steps);
s_manifold_pos_state = zeros(numsteps,size(state,2),manifold_steps);

us_manifold_neg_state = zeros(numsteps,size(state,2),manifold_steps);
s_manifold_neg_state = zeros(numsteps,size(state,2),manifold_steps);

us_manifold_pos_time = zeros(numsteps,manifold_steps);
us_manifold_neg_time= zeros(numsteps,manifold_steps);
s_manifold_pos_time= zeros(numsteps,manifold_steps);
s_manifold_neg_time= zeros(numsteps,manifold_steps);

us_manifold_pos_time_cross = zeros(10,manifold_steps);
us_manifold_pos_state_cross = zeros(10,size(state,2),manifold_steps);

us_manifold_neg_time_cross = zeros(10,manifold_steps);
us_manifold_neg_state_cross = zeros(10,size(state,2),manifold_steps);

s_manifold_pos_time_cross = zeros(10,manifold_steps);
s_manifold_pos_state_cross = zeros(10,size(state,2),manifold_steps);

s_manifold_neg_time_cross = zeros(10,manifold_steps);
s_manifold_neg_state_cross = zeros(10,size(state,2),manifold_steps);

bad_step = [];
% propogate the initial conditions to globalize the manifolds
% use an event function here to stop the propogation at the first poincare
% section crossing
tspan = linspace(0,manifold_tf,numsteps);
tspan_bw = linspace(-manifold_tf,0,numsteps);
for ii = 1:manifold_steps
    fprintf('Propogating manifold - %2g/%2g \n',ii,manifold_steps);
    try
     
    [t, state, t_cross,state_cross,~] = ...
        ode113(@(t,state)pcrtbp_ode(t,state,mu),tspan, us_manifold_pos_initial_state(ii,:),options_cross);
    us_manifold_pos_time(:,ii) = t;
    us_manifold_pos_state(:,:,ii) = state;
    us_manifold_pos_time_cross(1:length(t_cross),ii) = t_cross;
    us_manifold_pos_state_cross(1:length(t_cross),:,ii) = state_cross;
    
    [t, state, t_cross,state_cross,~] = ...
        ode113(@(t,state)pcrtbp_ode(t,state,mu),tspan, us_manifold_neg_initial_state(ii,:),options_cross);
    us_manifold_neg_time(:,ii) = t;
    us_manifold_neg_state(:,:,ii) = state;
    us_manifold_neg_time_cross(1:length(t_cross),ii) = t_cross;
    us_manifold_neg_state_cross(1:length(t_cross),:,ii) = state_cross;
    
    [t, state, t_cross,state_cross,~] = ...
        ode113(@(t,state)bw_pcrtbp_ode(t,state,mu),tspan_bw, s_manifold_pos_initial_state(ii,:),options_cross);
    s_manifold_pos_time(:,ii) = t;
    s_manifold_pos_state(:,:,ii) = state;
    s_manifold_pos_time_cross(1:length(t_cross),ii) = t_cross;
    s_manifold_pos_state_cross(1:length(t_cross),:,ii) = state_cross;
    
    [t, state, t_cross,state_cross,~] = ... 
        ode113(@(t,state)bw_pcrtbp_ode(t,state,mu),tspan_bw, s_manifold_neg_initial_state(ii,:),options_cross);
    s_manifold_neg_time(:,ii) = t;
    s_manifold_neg_state(:,:,ii) = state;
    s_manifold_neg_time_cross(1:length(t_cross),ii) = t_cross;
    s_manifold_neg_state_cross(1:length(t_cross),:,ii) = state_cross;
    
    catch err
        
        fprintf('Manifold propagation error. Try increasing error tolerance \n')
        bad_step = [bad_step ii];

        continue
    end
end

% remove any bad manifolds if there was an error
if ~isempty(bad_step)
    us_manifold_pos_state(:,:,bad_step) = [];
    us_manifold_neg_state(:,:,bad_step) = [];
    s_manifold_pos_state(:,:,bad_step) = [];
    s_manifold_neg_state(:,:,bad_step) = [];
    
    us_manifold_pos_time(:,bad_step) = [];
    us_manifold_neg_time(:,bad_step)= [];
    s_manifold_pos_time(:,bad_step)= [];
    s_manifold_neg_time(:,bad_step)= [];

    manifold_steps = size(us_manifold_pos_state,3);
end

manifold.uspm_initial = us_manifold_pos_initial_state;
manifold.spm_initial = s_manifold_pos_initial_state;
manifold.usnm_intial = us_manifold_neg_initial_state;
manifold.snm_initial = s_manifold_neg_initial_state;
manifold.periodic_initial = periodic_orbit_initial;

manifold.us_manifold_pos_time_cross = us_manifold_pos_time_cross;
manifold.us_manifold_neg_time_cross = us_manifold_neg_time_cross;
manifold.s_manifold_pos_time_cross = s_manifold_pos_time_cross;
manifold.s_manifold_neg_time_cross = s_manifold_neg_time_cross;

manifold.us_manifold_pos_state_cross = us_manifold_pos_state_cross;
manifold.us_manifold_neg_state_cross = us_manifold_neg_state_cross;
manifold.s_manifold_pos_state_cross = s_manifold_pos_state_cross;
manifold.s_manifold_neg_state_cross = s_manifold_neg_state_cross;

fprintf('Completed manifold globalization \n')
% plot manifolds if selected
switch constants.manifold_plot
    case 'true'
    fig_handle = figure;
    hold on;grid on; axis equal
    title('Invariant Manifolds')
    xlabel('X Axis')
    ylabel('Y Axis')
    
    plot(-mu,0,'k*')
    plot(1-mu,0,'k*')
    plot(state(:,1),state(:,2),'b')
    
    % plot the manifolds
    
    for ii = 1:manifold_steps
        uspms = us_manifold_pos_state(:,:,ii);
        usnms = us_manifold_neg_state(:,:,ii);
        
        spms = s_manifold_pos_state(:,:,ii);
        snms = s_manifold_neg_state(:,:,ii);
        
        plot(uspms(:,1),uspms(:,2),'r')
        plot(usnms(:,1),usnms(:,2),'r')
        plot(spms(:,1),spms(:,2),'g')
        plot(snms(:,1),snms(:,2),'g')
    end
    
     energy_contour(mu, per_orbit_E, fig_handle);
    case 'false'
        fprintf('No plotting \n')
end

function [value,isterminal,direction]=manifold_event(t,state,constants)

% check if manifold is crossing any of the poincare sections
mu = constants.mu;
section = constants.poincare_section;

% output the states
x = state(1);
y = state(2);

switch section
    case 1
       
        value = y; % -1 < x < 0
        
    case 2
         
           value = x-(1-mu); % y > 0

    case 3
        
        value = x-(1-mu); % y < 0
       
    case 4
        
        value = y; % x < -1
    case 'transfer'
        value = y; % x > 0
    case 'y_axis'
        value = x+mu;
        
end

isterminal = 0;
direction = 0;

