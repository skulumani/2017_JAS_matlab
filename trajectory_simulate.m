% function that will simulate a PCRTBP trajectory given an initial condition
% and time span
% should also plot the zero velocity energy curve

% september 4 2014
% added CRTBP STM propogation. need to add comment section at top

% september 8 2014
% modified backward propogation error - need to add negative to vector
% field
% 2 october 2014 - adding backward propogation option for CRTBP STM


function [t, state] = trajectory_simulate(x0, tspan,ode_model, constants)

mu = constants.mu;
% ode_model = constants.ode_model;
numsteps = constants.numsteps;

try
    ode_options = constants.ode_options;
catch err
    ode_options = odeset('RelTol',1e-9,'AbsTol',1e-9);     %set tolerences
end

tspan_fw = linspace(tspan(1), tspan(end),numsteps);
tspan_bw = linspace(tspan(1),tspan(2),numsteps);

% simulate the trajectory
switch ode_model
    case 'pcrtbp'
        [t,state]=ode113(@(t,state)pcrtbp_ode(t,state,mu),tspan_fw, x0,ode_options);
    case 'bw_pcrtbp'
        [t,state]=ode113(@(t,state)bw_pcrtbp_ode(t,state,mu),tspan_bw, x0,ode_options);
    case 'crtbp'
        [t,state]=ode113(@(t,state)crtbp_ode(t,state,mu),tspan_fw, x0,ode_options);
    case 'bw_crtbp'
        [t,state]=ode113(@(t,state)bw_crtbp_ode(t,state,mu),tspan_bw, x0,ode_options);
    case 'pcrtbp_stm'
        phi0 = eye(4,4);
        phi0=reshape(phi0,16,1);
        initial_stm = [x0;phi0];
        
        [t,state] = ode113(@(t,stm)pcrtbp_stm(t,stm,mu),tspan_fw,initial_stm,ode_options);
    case 'crtbp_stm'
        phi0 = eye(6,6);
        phi0=reshape(phi0,36,1);
        initial_stm = [x0;phi0];
        
        [t,state] = ode113(@(t,stm)crtbp_stm(t,stm,mu),tspan_fw,initial_stm,ode_options);
    case 'bw_crtbp_stm'
        phi0 = eye(6,6);
        phi0=reshape(phi0,36,1);
        initial_stm = [x0;phi0];
        
        [t,state] = ode113(@(t,stm)bw_crtbp_stm(t,stm,mu),tspan_bw,initial_stm,ode_options);

        
end


