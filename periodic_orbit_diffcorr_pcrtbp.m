function [x0_out, cross_time_out, phi_cross_out] = periodic_orbit_diffcorr_pcrtbp(initial_state,constants)
%PERIODIC_ORBIT_DIFFCORR Differential correction for periodic orbit
%
% [x0_out cross_time_out] = periodic_orbit_diffcorr(initial_state,constants)
%
% This will use differential correction to iteratively improve the initial
% condition of a periodic orbit in the PCRTBP. The goal is to modify the
% intial x position or initial y velocity to arrive at the half period
% location with no x velocity.
%
%   Inputs:
%       - initial_state - 4x1 initial state vector for nondimensional PCRTBP
%       - constants - structure variable with various integration
%       tolerances
%           - .RelTol - relative tolerance for Matlab ODE solver
%           - .AbsTol - absolute tolerance for Matlab ODE Solver
%           - .mu - mass paramater of system
%           - .tol - iteration tolerance for differential correction
%
%   Outputs:
%       - x0_out - 4x1 initial condition of the periodic orbit in the
%       rotating reference frame now corrected 
%       - cross_time_out - 1x1 nondimensional half period of the periodic
%       orbit
%       - phi_cross_out - 4x4 STM at half period of periodic orbit
%
%   Dependencies:
%       - libration_points.m - calculate location of lagrange points
%       - pcrtbp_eig.m - calculate the eigenvalues/eigenvectors for the
%       linearized system and determine the appropriate subspaces
%
%   Author:
%       - Shankar Kulumani 31 August 2014
%           - list revisions
%
% More detailed help is in the <a href="matlab: help
% periodic_orbit_initial>extended_help">extended help</a>.

% peform the differential correction to iterate to arrive at a better
% periodic orbit initial state

% integrate the nonlinear EOMS using the linear approximation intial state
x0 = initial_state;
TSPAN = [0 20] ;

RelTol = constants.RelTol;
AbsTol = constants.AbsTol;

tol = constants.tol;

mu = constants.mu;

diffcorr_plot=constants.diffcorr_plot;
diffcorr_select = constants.diffcorr_select_pcrtbp;

if diffcorr_plot == 1
    fig_handle = figure;
    hold all;grid on
    title('Differential Correction Iterations')
    xlabel('X Axis')
    ylabel('Y Axis')
end

options_cross = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events',@events_xcross);
options_stm = odeset('RelTol',RelTol,'AbsTol',AbsTol);

delta_xdot1 = 1;
iter = 1;
while abs(delta_xdot1) > tol && iter < 50
    fprintf('DiffCorr - %2g \n', iter)
    
    % high accuracy and Events ON for differential correction
    
    [~,~,cross_t,cross_state,ie] = ode113(@(t,state)pcrtbp_ode(t,state,mu),TSPAN,x0,options_cross) ;
    
    x1     = cross_state(end,1) ;
    y1     = cross_state(end,2) ;
    delta_xdot1  = cross_state(end,3) ;
    ydot1  = cross_state(end,4) ;
    
    % propogate the state transition matrix to the crossing time
    phi0 = eye(4,4);
    phi0=reshape(phi0,16,1);
    
    % augment the initial stm with the intial periodic orbit state for
    % propogation
    
    initial_stm = [x0;phi0];
    
    [t,stm] = ode113(@(t,stm)pcrtbp_stm(t,stm,mu),[0 cross_t(end)],initial_stm,options_stm);	% integration
    
    % extract state from 0 to crossing time
    state = stm(:,1:4);
    
    % pull out the final stm (phi(t1,t0) )
    phi1 = stm(end,5:end);
    phi1 = reshape(phi1,4,4);
    
    % calculate the state derivatives (f(x) at the crossing time)
    
    statedot=pcrtbp_ode(cross_t,cross_state,mu);
    
    xddot1 = statedot(3); % need the acceleration in x direction
    % compute the corrections (delta

    % select which initial state to modify to correct final crossing state (can
    % modify either x(0) or ydot(0)
    switch diffcorr_select
        case 'initial_ydot'
            delta_ydot = delta_xdot1 / ( phi1(3,4) - xddot1/ydot1*phi1(2,4));
            
            del = [0;0;0;delta_ydot];
        case 'initial_x'
            delta_x = delta_xdot1 / (phi1(3,1) - xddot1/ydot1*phi1(2,1));
            
            del = [delta_x;0;0;0];
    end
    
    x0 = x0-del;
    iter = iter+1;
    
    % switch to plot successive iterations on the intial condition
    
    if diffcorr_plot==1,
        set(0,'CurrentFigure',fig_handle)
        plot(state(:,1),state(:,2));
        m = length(state) ;
        plot(state(m,1),state(m,2),'bo');
        drawnow
    end
end % end of while loop to correct the initial condition

x0_out = x0;
cross_time_out = cross_t(end);
phi_cross_out = phi1; 

fprintf('DiffCorr Finished \n')
end

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
%
% Examples:
% foo(1,2,3)
%
% See also:
% BAR
% SOMECLASS/SOMEMETHOD
%
%   References
%       - Dynamical systems, the three-body problem and space mission
%       design Koon WS, Lo MW, Marsden JE, Ross SD.  World Scientific;
%       2000.
%

error('This is a placeholder function just for helptext');
end
