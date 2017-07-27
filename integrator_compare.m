% 27 Feb 2015 - script to compare performance of variational integrator and
% ODE45 and normal euler integration
function [] = integrator_compare()

addpath(genpath('./ode_solvers'));
constants = crtbp_constants;
x0 = [0.75;0;0.2;0]';
t0 = 0;
tf = 200; % figure out how to dimensionalize time
num_steps = 1 * 10 .^ [6, 7, 8];
width = 5;
height = 2;

% arrays to save the data
Ehist_ode45 = zeros(length(num_steps), 1e6);
Ehist_ode4 = zeros(length(num_steps), 1e6);
Ehist_trap = zeros(length(num_steps), 1e6);
meanE_ode45 = zeros(length(num_steps),1);
meanE_ode = zeros(length(num_steps), 1);
meanE_trap = zeros(length(num_steps),1);
cputime_ode45 = zeros(length(num_steps), 1);
cputime_ode = zeros(length(num_steps), 1);
cputime_trap = zeros(length(num_steps), 1);
step_size = zeros(length(num_steps), 1);

for ii = 1:length(num_steps)
    ns = num_steps(ii);
    t_vec = linspace(t0,tf,ns);
    h = (tf-t0)/ns;
    step_size(ii) = h;
    % propogate using ODE45
    ode45_st = tic;
    [~,state_ode45]=ode45(@(t,state)pcrtbp_ode(t,state,constants.mu),t_vec, x0, odeset('RelTol',1/ns,'AbsTol',1/ns));
    ode45_end = toc(ode45_st);
    
    % propogate using ODE4
    ode4_st = tic;
    state_ode4 =ode4(@(t,state)pcrtbp_ode(t,state,constants.mu),t_vec, x0);
    ode4_end = toc(ode4_st);
    
    % propogate using variational integrator rectangle
    % vir_st = tic;
    % [t_rect, state_rect] = pcrtbp_variational(x0,t0,tf, 'rect', ns,constants);
    % vir_end = toc(vir_st);
    
    % use trapezoidal mode
    vit_st = tic;
    [~, state_trap] = pcrtbp_variational(x0,t0,tf, 'trap', ns,constants);
    vit_end = toc(vit_st);
    
    % compute the energy difference
    % calculate jacobi energies for both
    [E_ode45] = energyconst(state_ode45(1:1000:end, :),constants.mu);
    [E_ode4] = energyconst(state_ode4(1:1000:end, :),constants.mu);
    % [E_rect] = energyconst(state_rect,constants.mu);
    [E_trap] = energyconst(state_trap(1:1000:end, :),constants.mu);
     
    % energy diff
    fprintf('\ntf = %5.2f nondim = %5.2e sec = %5.2f yrs\n', tf, tf*constants.t_scale, tf*constants.t_scale/86400/365);
    fprintf('N = %5.2e steps  \n', ns);
    fprintf('h = %5.2e nondim = %5.2f sec = %5.2f days\n\n',h,h*constants.t_scale, h*constants.t_scale/86400);
    
    fprintf('ODE45 energy drift %12.10e\n',mean(abs(E_ode45(end)-E_ode45(1))));
    fprintf('ODE4 energy drift %12.10e\n',mean(abs(E_ode4(end)-E_ode4(1))));
    % fprintf('VI RECT energy drift %12.10e\n',mean(abs(E_rect-E_rect(1))));
    fprintf('VI TRAP Energy drift %12.10e\n',mean(abs(E_trap-E_trap(1))));
    
    fprintf('ODE45 mean E var %12.10e\n', var(E_ode45));
    fprintf('ODE4 mean E var %12.10e\n', var(E_ode4));
    % fprintf('VI RECT mean E var %12.10e\n',var(E_rect));
    fprintf('VI TRAP mean E var %12.10e\n',var(E_trap));
    
    fprintf('ODE45 run time %5.2f sec\n',ode45_end);
    fprintf('ODE4 run time %5.2f sec\n',ode4_end);
    % fprintf('VI RECT run time %5.2f sec\n',vir_end);
    fprintf('VI TRAP run time %5.2f sec\n',vit_end);
    
    % store in arrays
    Ehist_ode45(ii,1:length(E_ode45)) = E_ode45;
    Ehist_ode4(ii, 1:length(E_ode4)) = E_ode4;
    Ehist_trap(ii,1:length(E_trap)) = E_trap;

    meanE_ode45(ii) = mean(abs(E_ode45 - E_ode45(1)));
    meanE_ode4(ii) = mean(abs(E_ode4 - E_ode4(1)));
    meanE_trap(ii) = mean(abs(E_trap - E_trap(1)));
    
    clear E_ode45 E_ode4 E_trap state_trap state_ode45 state_ode4

    cputime_ode45(ii) = ode45_end;
    cputime_ode4(ii) = ode4_end;
    cputime_trap(ii) = vit_end;
end

save('E_comparison.mat', 'Ehist_ode45', 'Ehist_ode4', 'Ehist_trap', 'meanE_ode45', 'meanE_ode4', 'meanE_trap', 'cputime_ode45', 'cputime_ode4', 'cputime_trap', 'num_steps', 'step_size', '-v7.3')
